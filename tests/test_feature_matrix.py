from pathlib import Path

import pytest
from pandas import DataFrame

from flash_fraq.feature_matrix import (
    ABLATION_ONLY_COLUMNS,
    REGRESSION_LABEL_COLUMN,
    FeatureMatrixAssembler,
    FeatureMatrixConfig,
    pivot_cnv_bins,
)

SAMPLE_IDS = ["mix_001", "mix_002"]


def _write_sample_metadata(input_dir: Path) -> None:
    DataFrame(
        {
            "sample_id": SAMPLE_IDS,
            "tumor_source_id": ["HCC1395_a", "HCC1395_a"],
            "healthy_source_id": ["BH01", "IA01"],
            "platform_tumor": ["NovaSeq", "NovaSeq"],
            "platform_healthy": ["HiSeq2500", "HiSeq2500"],
            "reference_build": ["hg19", "hg19"],
            "target_coverage": [0.5, 1.0],
            "nominal_tumor_fraction": [0.10, 0.20],
            "empirical_tumor_fraction": [0.083, 0.174],
            "tumor_subject_id": ["HCC1395", "HCC1395"],
            "healthy_subject_id": ["BH01_subject", "IA01_subject"],
            "mtdna_fraction_log10": [-2.1, -2.4],
        }
    ).to_csv(input_dir / "sample_metadata.csv", index=False)


def _write_fragment_stats(input_dir: Path) -> None:
    DataFrame(
        {
            "sample_id": SAMPLE_IDS,
            "fragment_size_median": [166.0, 168.0],
            "fragment_size_mode": [167.0, 167.0],
            "fragment_size_entropy": [4.21, 4.19],
        }
    ).to_csv(input_dir / "fragment_stats.csv", index=False)


ARM_LABELS = ("chr1p", "chr1q", "chr2p")


def _write_arm_zscores(input_dir: Path) -> None:
    rows = [
        {"sample_id": sample_id, "arm": arm, "zscore": 0.1 * (index + 1)}
        for sample_id in SAMPLE_IDS
        for index, arm in enumerate(ARM_LABELS)
    ]
    DataFrame(rows).to_csv(input_dir / "arm_zscores.csv", index=False)


def _write_window_ratios(input_dir: Path) -> None:
    rows = [
        {
            "sample_id": sample_id,
            "window_index": window_index,
            "short_long_ratio": 0.5 + 0.01 * window_index,
            "gc_corrected_total_coverage": 10.0 + window_index,
            "gc_corrected_short_coverage": 4.0 + window_index,
        }
        for sample_id in SAMPLE_IDS
        for window_index in range(1, 4)
    ]
    DataFrame(rows).to_csv(input_dir / "window_ratios.csv", index=False)


def _write_cnv_bins(input_dir: Path) -> None:
    rows = [
        {"sample_id": sample_id, "bin_index": bin_index, "cnv_logr": 0.02 * bin_index}
        for sample_id in SAMPLE_IDS
        for bin_index in range(1, 3)
    ]
    DataFrame(rows).to_csv(input_dir / "cnv_bins.csv", index=False)


def _write_ichor_estimates(input_dir: Path) -> None:
    DataFrame(
        {
            "sample_id": SAMPLE_IDS,
            "ichorCNA_tumor_fraction_estimate": [0.09, 0.19],
        }
    ).to_csv(input_dir / "ichor_estimates.csv", index=False)


def _write_split_assignment(input_dir: Path) -> None:
    DataFrame(
        {
            "subject_id": ["HCC1395", "BH01_subject", "IA01_subject"],
            "split_group": ["train", "train", "train"],
        }
    ).to_csv(input_dir / "split_assignment.csv", index=False)


def _build_assembler(tmp_path: Path) -> FeatureMatrixAssembler:
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    _write_sample_metadata(input_dir)
    _write_fragment_stats(input_dir)
    _write_arm_zscores(input_dir)
    _write_window_ratios(input_dir)
    _write_cnv_bins(input_dir)
    _write_ichor_estimates(input_dir)
    _write_split_assignment(input_dir)

    config = FeatureMatrixConfig.from_sources(
        input_dir=input_dir,
        output_path=tmp_path / "outputs" / "feature_matrix.csv",
        # This fixture's cnv_bins.csv (_write_cnv_bins) carries no chrom/start -- it
        # predates the centromere/chromosome-end mask and isn't testing it. Masking
        # itself is covered by test_pivot_cnv_bins_masks_bins_near_centromere_or_end
        # below and tests/test_genome_mask.py.
        mask_margin_bp=0,
    )
    return FeatureMatrixAssembler(config)


def test_assemble_produces_one_row_per_sample(tmp_path: Path) -> None:
    assembler = _build_assembler(tmp_path)

    matrix = assembler.assemble()

    assert matrix.shape[0] == len(SAMPLE_IDS)
    assert sorted(matrix["sample_id"]) == sorted(SAMPLE_IDS)


def test_assemble_includes_all_schema_columns(tmp_path: Path) -> None:
    assembler = _build_assembler(tmp_path)

    matrix = assembler.assemble()

    expected_columns = {
        "sample_id",
        "tumor_source_id",
        "healthy_source_id",
        "platform_tumor",
        "platform_healthy",
        "reference_build",
        "target_coverage",
        "nominal_tumor_fraction",
        "empirical_tumor_fraction",
        "fragment_size_median",
        "fragment_size_mode",
        "fragment_size_entropy",
        "ichorCNA_tumor_fraction_estimate",
        "mtdna_fraction_log10",
        "split_group",
    }
    expected_columns |= {f"short_long_ratio_window_{i}" for i in range(1, 4)}
    expected_columns |= {f"gc_corrected_total_coverage_window_{i}" for i in range(1, 4)}
    expected_columns |= {f"gc_corrected_short_coverage_window_{i}" for i in range(1, 4)}
    expected_columns |= {f"arm_zscore_{arm}" for arm in ARM_LABELS}
    expected_columns |= {f"cnv_logr_bin_{i}" for i in range(1, 3)}

    assert expected_columns.issubset(set(matrix.columns))


def test_empirical_tumor_fraction_is_the_regression_label(tmp_path: Path) -> None:
    assembler = _build_assembler(tmp_path)

    matrix = assembler.assemble()

    assert REGRESSION_LABEL_COLUMN == "empirical_tumor_fraction"
    assert "empirical_tumor_fraction" in matrix.columns
    assert matrix.loc[matrix["sample_id"] == "mix_001", "empirical_tumor_fraction"].iloc[0] == 0.083


def test_ichorcna_estimate_is_flagged_ablation_only_and_excluded_from_primary_features(
    tmp_path: Path,
) -> None:
    assembler = _build_assembler(tmp_path)

    matrix = assembler.assemble()
    primary_features = assembler.get_primary_feature_columns(matrix)

    assert "ichorCNA_tumor_fraction_estimate" in ABLATION_ONLY_COLUMNS
    assert "ichorCNA_tumor_fraction_estimate" in matrix.columns
    assert "ichorCNA_tumor_fraction_estimate" not in primary_features


def test_mtdna_fraction_log10_stays_metadata_excluded_from_primary_features(tmp_path: Path) -> None:
    """audit finding 2: mtdna_fraction_log10 is a tumor-read counter (HCC1395 gDNA
    ~10^-2.6 vs cfDNA ~10^-4.9), not a fragmentomics feature -- it stays in the
    matrix as metadata but is never a candidate feature."""
    assembler = _build_assembler(tmp_path)

    matrix = assembler.assemble()
    primary_features = assembler.get_primary_feature_columns(matrix)

    assert "mtdna_fraction_log10" in matrix.columns
    assert "mtdna_fraction_log10" not in primary_features


def test_split_group_assigned_at_source_subject_level_before_mixing(tmp_path: Path) -> None:
    assembler = _build_assembler(tmp_path)

    matrix = assembler.assemble()

    # Both mixed samples draw from source subjects assigned to the "train" split_group
    # (assigned at the subject level, before mixing) so both rows inherit "train".
    assert set(matrix["split_group"]) == {"train"}


def test_mismatched_source_subject_splits_raise(tmp_path: Path) -> None:
    assembler = _build_assembler(tmp_path)

    sample_metadata = assembler.load_sample_metadata()
    mismatched_split_assignment = DataFrame(
        {
            "subject_id": ["HCC1395", "BH01_subject", "IA01_subject"],
            "split_group": ["train", "train", "val"],
        }
    )

    try:
        assembler.assign_split_groups(sample_metadata, mismatched_split_assignment)
    except ValueError as error:
        assert "mix_002" in str(error)
    else:
        raise AssertionError("expected ValueError for mismatched source-subject split groups")


def test_write_creates_output_file(tmp_path: Path) -> None:
    assembler = _build_assembler(tmp_path)
    matrix = assembler.assemble()

    output_path = assembler.write(matrix)

    assert output_path.exists()
    assert output_path == assembler.config.output_path


def test_pivot_cnv_bins_masks_bins_near_centromere_or_end() -> None:
    """docs/results-corruption-audit.md's coverage-bias finding: cnv_logr_bin_* columns
    for a bin within the default 5 Mb margin of chr21's centromere (11,288,129-
    14,288,129, hg19_centromeres.tsv) or either chromosome end (48,129,895 bp,
    hg19_chrom_sizes.tsv) must not appear in the pivoted matrix."""
    long_df = DataFrame(
        {
            "sample_id": ["mix_001"] * 4,
            "bin_index": [0, 1, 2, 3],
            "chrom": ["chr21"] * 4,
            # 0: chromosome-start-adjacent; 1: centromere-adjacent; 2: chromosome-end-
            # adjacent; 3: mid chr21q, far from all three -- the one survivor.
            "start": [0, 12_000_000, 47_000_000, 25_000_000],
            "cnv_logr": [0.1, 0.2, 0.3, 0.4],
        }
    )

    wide = pivot_cnv_bins(long_df)

    assert list(wide.columns) == ["cnv_logr_bin_3"]
    assert wide["cnv_logr_bin_3"].iloc[0] == pytest.approx(0.4)


def test_pivot_cnv_bins_mask_margin_bp_zero_disables_masking() -> None:
    long_df = DataFrame(
        {
            "sample_id": ["mix_001"] * 2,
            "bin_index": [0, 1],
            "chrom": ["chr21"] * 2,
            "start": [0, 25_000_000],
            "cnv_logr": [0.1, 0.4],
        }
    )

    wide = pivot_cnv_bins(long_df, mask_margin_bp=0)

    assert sorted(wide.columns) == ["cnv_logr_bin_0", "cnv_logr_bin_1"]
