import json
from pathlib import Path

import numpy as np
import pytest
from pandas import DataFrame, concat

from flash_fraq.decoy_qc import DecoyClassifierQc, DecoyClassifierQcConfig

N_WINDOWS = 5
SOURCES = ("GSE71378", "Cristiano", "HCC1395")


def _build_feature_matrix(
    n_per_source: int,
    leak_platform_into_fraction: bool,
    random_state: int,
    regime: str = "real_background",
    library_confound_rho: float = 0.0,
) -> DataFrame:
    """Synthetic feature matrix with an injected, easily-separable platform signal.

    If `leak_platform_into_fraction` is True, tumor fraction is assigned so that it
    correlates with source of origin (the batch-effect failure mode §6 warns about).
    Otherwise tumor fraction is assigned independently of source.

    `library_confound_rho` injects a non-CNV feature that tracks `nominal_tumor_fraction`
    directly (check 3's failure mode) independent of `leak_platform_into_fraction`
    (check 2's), so the two checks can be exercised independently.
    """
    rng = np.random.default_rng(random_state)
    rows = []
    sample_index = 0
    for source_offset, source in enumerate(SOURCES):
        for _ in range(n_per_source):
            if leak_platform_into_fraction:
                tumor_fraction = float(np.clip(rng.normal(loc=source_offset * 0.3, scale=0.05), 0.0, 1.0))
            else:
                tumor_fraction = float(rng.uniform(0.0, 1.0))

            row = {
                "sample_id": f"sample_{sample_index}",
                "tumor_source_id": "HCC1395",
                "healthy_source_id": source if source != "HCC1395" else "GSE71378",
                "platform_tumor": "centerA",
                "platform_healthy": source,
                "empirical_tumor_fraction": tumor_fraction,
                "nominal_tumor_fraction": tumor_fraction,
                "ichorCNA_tumor_fraction_estimate": tumor_fraction + rng.normal(scale=0.02),
                "split_group": "train",
                "regime": regime,
            }
            for window in range(1, N_WINDOWS + 1):
                # Platform-specific offset injected into every feature column, independent
                # of tumor fraction unless leak_platform_into_fraction ties source to fraction.
                row[f"short_long_ratio_window_{window}"] = source_offset * 2.0 + rng.normal(scale=0.1)
                row[f"gc_corrected_total_coverage_window_{window}"] = source_offset * 1.5 + rng.normal(scale=0.1)
                row[f"arm_zscore_{window}"] = source_offset * 1.0 + rng.normal(scale=0.1)
                row[f"cnv_logr_bin_{window}"] = source_offset * 0.5 + rng.normal(scale=0.1)
            row["fragment_size_median"] = source_offset * 3.0 + rng.normal(scale=0.1)
            row["fragment_size_mode"] = source_offset * 3.0 + rng.normal(scale=0.1)
            row["fragment_size_entropy"] = source_offset * 0.2 + rng.normal(scale=0.1)
            row["mtdna_fraction_log10"] = rng.normal(scale=0.1)
            if library_confound_rho:
                row["fragment_size_short_ratio_library"] = library_confound_rho * tumor_fraction + rng.normal(
                    scale=0.02
                )

            rows.append(row)
            sample_index += 1

    return DataFrame(rows)


def _run_qc(tmp_path: Path, leak_platform_into_fraction: bool) -> tuple[bool, float]:
    feature_matrix = _build_feature_matrix(
        n_per_source=300, leak_platform_into_fraction=leak_platform_into_fraction, random_state=0
    )
    feature_matrix_path = tmp_path / "feature_matrix.csv"
    feature_matrix.to_csv(feature_matrix_path, index=False)

    config = DecoyClassifierQcConfig.from_sources(
        feature_matrix_path=feature_matrix_path,
        output_dir=tmp_path / "qc_output",
    )
    report = DecoyClassifierQc(config).run()
    return report.passed, report.correlation_check.max_abs_correlation


def test_decoy_qc_flags_platform_correlated_with_tumor_fraction(tmp_path: Path) -> None:
    passed, max_abs_correlation = _run_qc(tmp_path, leak_platform_into_fraction=True)

    assert passed is False
    assert max_abs_correlation > 0.1
    assert (tmp_path / "qc_output" / "decoy_qc_report.json").exists()


def test_decoy_qc_passes_when_tumor_fraction_independent_of_platform(tmp_path: Path) -> None:
    passed, max_abs_correlation = _run_qc(tmp_path, leak_platform_into_fraction=False)

    assert passed is True
    assert max_abs_correlation <= 0.1
    assert (tmp_path / "qc_output" / "decoy_qc_report.json").exists()


def test_regime_scoped_runs_do_not_overwrite_each_others_report(tmp_path: Path) -> None:
    # Both admixture regimes are commonly run into the same --output-dir (one process
    # per --regime); the report filename must be regime-scoped or the second run
    # silently clobbers the first's report.
    output_dir = tmp_path / "qc_output"
    reports = {}
    for regime in ("real_background", "matched_pair"):
        feature_matrix = _build_feature_matrix(
            n_per_source=30, leak_platform_into_fraction=False, random_state=2, regime=regime
        )
        feature_matrix_path = tmp_path / f"feature_matrix_{regime}.csv"
        feature_matrix.to_csv(feature_matrix_path, index=False)

        config = DecoyClassifierQcConfig(
            feature_matrix_path=feature_matrix_path,
            output_dir=output_dir,
            regime=regime,
        )
        qc = DecoyClassifierQc(config)
        qc.run()
        reports[regime] = qc.report_path()

    assert reports["real_background"] != reports["matched_pair"]
    assert reports["real_background"].exists()
    assert reports["matched_pair"].exists()


def test_fixed_fraction_check_detects_injected_platform_signal(tmp_path: Path) -> None:
    feature_matrix = _build_feature_matrix(n_per_source=60, leak_platform_into_fraction=False, random_state=1)
    feature_matrix_path = tmp_path / "feature_matrix.csv"
    feature_matrix.to_csv(feature_matrix_path, index=False)

    config = DecoyClassifierQcConfig(
        feature_matrix_path=feature_matrix_path,
        output_dir=tmp_path / "qc_output",
        fixed_fraction_target=0.5,
        fixed_fraction_tolerance=1.0,
    )
    result = DecoyClassifierQc(config).run_fixed_fraction_check(feature_matrix)

    assert result.source_discriminable is True
    assert result.accuracy > result.baseline_accuracy


COHORT_SOURCES = (("Snyder", "HiSeq2000"), ("Cristiano", "HiSeq2500"))


def _build_cohort_matrix(n_per_source: int, leak_platform_into_fraction: bool, random_state: int) -> DataFrame:
    """Single-source (pre-admixture) matrix, as cohort_matrix.py assembles it.

    Also carries the two things a matrix built from real per-sample tables carries and
    the synthetic §9 matrix above does not: a feature column missing for one whole
    source, and a constant column.
    """
    rng = np.random.default_rng(random_state)
    rows = []
    for source_offset, (dataset, platform) in enumerate(COHORT_SOURCES):
        for sample_index in range(n_per_source):
            if leak_platform_into_fraction:
                tumor_fraction = float(np.clip(rng.normal(loc=source_offset * 0.3, scale=0.05), 0.0, 1.0))
            else:
                tumor_fraction = float(rng.uniform(0.0, 1.0))

            row = {
                "sample_id": f"{dataset}_{sample_index}",
                "source_dataset": dataset,
                "platform": platform,
                "reference_build": "hg19",
                "ichorCNA_tumor_fraction_estimate": tumor_fraction,
                "fragment_size_median": source_offset * 3.0 + rng.normal(scale=0.1),
                # Present only for the second source: must be dropped, not imputed into a
                # perfect source giveaway.
                "short_long_ratio_window_99": rng.normal() if source_offset else float("nan"),
                "cnv_logr_bin_99": 1.0,
            }
            for window in range(1, N_WINDOWS + 1):
                row[f"short_long_ratio_window_{window}"] = source_offset * 2.0 + rng.normal(scale=0.1)
                row[f"arm_zscore_chr{window}p"] = source_offset * 1.0 + rng.normal(scale=0.1)
                # Two orders of magnitude above every other column, as the real
                # GC-corrected coverages are.
                row[f"gc_corrected_total_coverage_window_{window}"] = 50_000.0 + rng.normal(scale=500.0)
            rows.append(row)

    return DataFrame(rows)


def _run_cohort_qc(tmp_path: Path, leak_platform_into_fraction: bool):
    feature_matrix = _build_cohort_matrix(
        n_per_source=200, leak_platform_into_fraction=leak_platform_into_fraction, random_state=0
    )
    feature_matrix_path = tmp_path / "cohort_feature_matrix.csv"
    feature_matrix.to_csv(feature_matrix_path, index=False)

    config = DecoyClassifierQcConfig.for_cohort(
        feature_matrix_path=feature_matrix_path, output_dir=tmp_path / "qc_output"
    )
    return DecoyClassifierQc(config).run()


def test_cohort_qc_flags_source_correlated_with_ichor_fraction(tmp_path: Path) -> None:
    report = _run_cohort_qc(tmp_path, leak_platform_into_fraction=True)

    assert report.passed is False
    assert report.correlation_check.max_abs_correlation > 0.1
    assert set(report.correlation_check.correlation_by_source) == {"Snyder|HiSeq2000", "Cristiano|HiSeq2500"}


def test_cohort_qc_passes_when_fraction_independent_of_source(tmp_path: Path) -> None:
    report = _run_cohort_qc(tmp_path, leak_platform_into_fraction=False)

    assert report.passed is True
    assert report.correlation_check.max_abs_correlation <= 0.1


def test_correlation_check_records_its_folds_and_the_direct_source_correlation(tmp_path: Path) -> None:
    """The out-of-fold number must be reproducible from the report alone, and readable
    against the plain corr(source indicator, fraction) it could collapse onto."""
    report = _run_cohort_qc(tmp_path, leak_platform_into_fraction=True)
    check = report.correlation_check

    assert check.skipped_reason is None
    assert check.n_splits == 5
    assert check.fold_random_state == 0
    assert set(check.source_indicator_correlation_by_source) == set(check.correlation_by_source)
    assert check.max_abs_source_indicator_correlation > 0.1


def test_correlation_check_skips_when_a_source_is_smaller_than_n_splits(tmp_path: Path) -> None:
    """Fewer samples in a source than folds is a precondition failure, not a clean pass:
    sklearn would raise, and reporting 0.0 would read as measured evidence of no leakage."""
    feature_matrix = _build_cohort_matrix(n_per_source=20, leak_platform_into_fraction=True, random_state=4)
    # Snyder drops to 3 samples, below the 5 folds check 2 asks for.
    snyder = feature_matrix.index[feature_matrix["source_dataset"] == "Snyder"]
    feature_matrix = feature_matrix.drop(index=snyder[3:])
    feature_matrix_path = tmp_path / "cohort_feature_matrix.csv"
    feature_matrix.to_csv(feature_matrix_path, index=False)
    config = DecoyClassifierQcConfig.for_cohort(
        feature_matrix_path=feature_matrix_path, output_dir=tmp_path / "qc_output"
    )

    report = DecoyClassifierQc(config).run()
    check = report.correlation_check

    assert check.max_abs_correlation is None
    assert check.flagged is False
    assert check.n_splits is None
    assert "3 samples" in check.skipped_reason
    assert "n_splits=5" in check.skipped_reason
    # A skipped check writes null, not a number a reader could take for a measurement.
    written = json.loads((tmp_path / "qc_output" / "decoy_qc_report.json").read_text())
    assert written["correlation_check"]["max_abs_correlation"] is None


def test_correlation_check_runs_with_fewer_folds_when_configured(tmp_path: Path) -> None:
    """The guard is on n_splits, not on a hardcoded 5: a smaller fold count runs."""
    feature_matrix = _build_cohort_matrix(n_per_source=20, leak_platform_into_fraction=True, random_state=4)
    snyder = feature_matrix.index[feature_matrix["source_dataset"] == "Snyder"]
    feature_matrix = feature_matrix.drop(index=snyder[3:])
    feature_matrix_path = tmp_path / "cohort_feature_matrix.csv"
    feature_matrix.to_csv(feature_matrix_path, index=False)
    config = DecoyClassifierQcConfig.for_cohort(
        feature_matrix_path=feature_matrix_path, output_dir=tmp_path / "qc_output", n_splits=3
    )

    check = DecoyClassifierQc(config).run().correlation_check

    assert check.skipped_reason is None
    assert check.n_splits == 3
    assert check.max_abs_correlation is not None


def test_ragged_and_constant_feature_columns_are_dropped(tmp_path: Path) -> None:
    feature_matrix = _build_cohort_matrix(n_per_source=20, leak_platform_into_fraction=False, random_state=2)
    feature_matrix_path = tmp_path / "cohort_feature_matrix.csv"
    feature_matrix.to_csv(feature_matrix_path, index=False)
    config = DecoyClassifierQcConfig.for_cohort(
        feature_matrix_path=feature_matrix_path, output_dir=tmp_path / "qc_output"
    )

    features = DecoyClassifierQc(config).prepare_features(feature_matrix)

    assert "short_long_ratio_window_99" not in features.columns
    assert "cnv_logr_bin_99" not in features.columns
    assert "short_long_ratio_window_1" in features.columns
    assert not features.isna().to_numpy().any()
    # Metadata and the stand-in label never reach the classifier.
    assert {"sample_id", "source_dataset", "platform", "reference_build"}.isdisjoint(features.columns)
    assert "ichorCNA_tumor_fraction_estimate" not in features.columns


def test_library_confound_check_flags_a_non_cnv_feature_tracking_nominal_fraction(tmp_path: Path) -> None:
    feature_matrix = _build_feature_matrix(
        n_per_source=300, leak_platform_into_fraction=False, random_state=5, library_confound_rho=3.0
    )
    feature_matrix_path = tmp_path / "feature_matrix.csv"
    feature_matrix.to_csv(feature_matrix_path, index=False)
    config = DecoyClassifierQcConfig.from_sources(feature_matrix_path=feature_matrix_path, output_dir=tmp_path / "qc")

    report = DecoyClassifierQc(config).run()
    check = report.library_confound_check

    assert check.skipped_reason is None
    assert check.spearman_rho is not None
    assert abs(check.spearman_rho) > 0.1
    assert check.flagged is True
    assert report.passed is False
    assert abs(check.per_feature_spearman_rho["fragment_size_short_ratio_library"]) > 0.1


def test_library_confound_check_passes_when_fragment_features_are_independent_of_fraction(tmp_path: Path) -> None:
    feature_matrix = _build_feature_matrix(n_per_source=300, leak_platform_into_fraction=False, random_state=0)
    feature_matrix_path = tmp_path / "feature_matrix.csv"
    feature_matrix.to_csv(feature_matrix_path, index=False)
    config = DecoyClassifierQcConfig.from_sources(feature_matrix_path=feature_matrix_path, output_dir=tmp_path / "qc")

    report = DecoyClassifierQc(config).run()
    check = report.library_confound_check

    assert check.skipped_reason is None
    assert check.residualised is True
    assert abs(check.spearman_rho) <= 0.1
    assert check.flagged is False
    assert report.passed is True


def test_pooled_regime_matrix_refuses_to_run_without_an_explicit_regime(tmp_path: Path) -> None:
    real_background = _build_feature_matrix(
        n_per_source=20, leak_platform_into_fraction=False, random_state=7, regime="real_background"
    )
    matched_pair = _build_feature_matrix(
        n_per_source=5, leak_platform_into_fraction=False, random_state=8, regime="matched_pair"
    )
    pooled = concat([real_background, matched_pair], ignore_index=True)
    feature_matrix_path = tmp_path / "feature_matrix.csv"
    pooled.to_csv(feature_matrix_path, index=False)
    config = DecoyClassifierQcConfig(feature_matrix_path=feature_matrix_path, output_dir=tmp_path / "qc")

    with pytest.raises(ValueError, match="pools 2 regimes"):
        DecoyClassifierQc(config).run()


def test_explicit_regime_filters_the_matrix_before_running(tmp_path: Path) -> None:
    real_background = _build_feature_matrix(
        n_per_source=100, leak_platform_into_fraction=False, random_state=7, regime="real_background"
    )
    matched_pair = _build_feature_matrix(
        n_per_source=5, leak_platform_into_fraction=False, random_state=8, regime="matched_pair"
    )
    pooled = concat([real_background, matched_pair], ignore_index=True)
    feature_matrix_path = tmp_path / "feature_matrix.csv"
    pooled.to_csv(feature_matrix_path, index=False)
    config = DecoyClassifierQcConfig(
        feature_matrix_path=feature_matrix_path, output_dir=tmp_path / "qc", regime="real_background"
    )

    report = DecoyClassifierQc(config).run()

    assert report.library_confound_check.n_samples == len(real_background)


def test_constant_label_column_refuses_to_report_a_pass(tmp_path: Path) -> None:
    feature_matrix = _build_cohort_matrix(n_per_source=20, leak_platform_into_fraction=False, random_state=3)
    feature_matrix["ichorCNA_tumor_fraction_estimate"] = 0.0
    feature_matrix_path = tmp_path / "cohort_feature_matrix.csv"
    feature_matrix.to_csv(feature_matrix_path, index=False)
    config = DecoyClassifierQcConfig.for_cohort(
        feature_matrix_path=feature_matrix_path, output_dir=tmp_path / "qc_output"
    )

    with pytest.raises(ValueError, match="no variance"):
        DecoyClassifierQc(config).run()
