"""Cohort feature matrix (pre-admixture Snyder + Cristiano join).

The per-sample tables are produced by driving the *real* GcCorrector rather than by
hand-writing arm_zscores.csv, so the cohort pivot is tested against the arm labels
gc_correction actually emits, not against a restatement of them.
"""

from pathlib import Path

import numpy as np
import pytest
from pandas import DataFrame, read_csv

from flash_fraq.cohort_matrix import (
    COHORT_TUMOR_FRACTION_COLUMN,
    CohortFeatureMatrixAssembler,
    CohortMatrixConfig,
)
from flash_fraq.gc_correction import GcCorrectionConfig, GcCorrector

WINDOW_BP = 5_000_000
# chr1's centromere sits at 123.4 Mb, chr2's at 93.9 Mb, so 40 windows of chr1 and 30 of
# chr2 span all four of chr1p/chr1q/chr2p/chr2q.
CHROM_WINDOW_COUNTS = {"chr1": 40, "chr2": 30}


def _window_counts(sample_id: str, chrom_window_counts: dict[str, int], random_state: int) -> DataFrame:
    rng = np.random.default_rng(random_state)
    rows = []
    window_index = 1
    for chrom, n_windows in chrom_window_counts.items():
        for window in range(n_windows):
            start = window * WINDOW_BP
            rows.append(
                {
                    "sample_id": sample_id,
                    # Stable across samples: it indexes the shared window BED, which is
                    # what makes a cohort-wide pivot on window_index meaningful.
                    "window_index": window_index,
                    "chrom": chrom,
                    "start": start,
                    "end": start + WINDOW_BP,
                    "gc_content": float(rng.uniform(0.35, 0.55)),
                    "short_count": int(rng.integers(5_000, 15_000)),
                    "long_count": int(rng.integers(5_000, 15_000)),
                }
            )
            window_index += 1
    return DataFrame(rows)


def _write_sample(
    derived_root: Path,
    dataset: str,
    sample_id: str,
    chrom_window_counts: dict[str, int],
    tumor_fraction: float,
    random_state: int,
    complete: bool = True,
) -> None:
    sample_dir = derived_root / dataset / sample_id
    (sample_dir / "ichor").mkdir(parents=True, exist_ok=True)

    counts_path = sample_dir / "window_fragment_counts.tsv"
    _window_counts(sample_id, chrom_window_counts, random_state).to_csv(counts_path, sep="\t", index=False)
    GcCorrector(
        GcCorrectionConfig.from_sources(
            input_windows_path=counts_path,
            output_windows_path=sample_dir / "window_fragment_counts.gc_corrected.tsv",
            # This fixture's windows sit at real chr1/chr2 starts (window_index 1 and 41
            # are start=0) purely to get realistic arm buckets from HG19_CENTROMERE_MB --
            # they aren't testing the centromere/chromosome-end mask, which would drop
            # those first windows and shift the fixed window_index columns this module's
            # tests assert on. Masking itself is covered by tests/test_gc_correction.py
            # and tests/test_genome_mask.py.
            mask_margin_bp=0,
        )
    ).run_full()

    DataFrame(
        [{"sample_id": sample_id, "fragment_size_median": 166.0, "fragment_size_mode": 167, "fragment_size_entropy": 7.5}]
    ).to_csv(sample_dir / "fragment_stats.csv", index=False)

    DataFrame(
        [
            {"sample_id": sample_id, "bin_index": bin_index, "chrom": "chr1", "start": bin_index * 1_000_000, "cnv_logr": 0.1 * bin_index}
            for bin_index in range(1, 6)
        ]
    ).to_csv(sample_dir / "ichor" / "cnv_bins.csv", index=False)

    DataFrame([{"sample_id": sample_id, COHORT_TUMOR_FRACTION_COLUMN: tumor_fraction}]).to_csv(
        sample_dir / "ichor" / "ichor_estimates.csv", index=False
    )

    if complete:
        (sample_dir / "localfeat_sample_complete.flag").touch()


def _build_cohort(derived_root: Path) -> None:
    """Two datasets; SNY02 is missing all of chr2, so it has two arms fewer than the rest."""
    _write_sample(derived_root, "Snyder", "SNY01", CHROM_WINDOW_COUNTS, tumor_fraction=0.0, random_state=0)
    _write_sample(derived_root, "Snyder", "SNY02", {"chr1": 40}, tumor_fraction=0.1, random_state=1)
    _write_sample(derived_root, "Cristiano", "CRI01", CHROM_WINDOW_COUNTS, tumor_fraction=0.2, random_state=2)


def _config(derived_root: Path, output_path: Path) -> CohortMatrixConfig:
    return CohortMatrixConfig.from_sources(
        derived_root=derived_root,
        dataset_platforms={"Snyder": "HiSeq2000", "Cristiano": "HiSeq2500"},
        output_path=output_path,
    )


def test_assembles_one_row_per_sample_tagged_by_source(tmp_path: Path) -> None:
    _build_cohort(tmp_path / "derived")
    assembler = CohortFeatureMatrixAssembler(_config(tmp_path / "derived", tmp_path / "cohort.csv"))

    matrix = assembler.assemble()

    assert sorted(matrix["sample_id"]) == ["CRI01", "SNY01", "SNY02"]
    assert dict(zip(matrix["sample_id"], matrix["source_dataset"])) == {
        "SNY01": "Snyder",
        "SNY02": "Snyder",
        "CRI01": "Cristiano",
    }
    assert set(matrix["platform"]) == {"HiSeq2000", "HiSeq2500"}
    assert set(matrix["reference_build"]) == {"hg19"}
    assert matrix.loc[matrix["sample_id"] == "CRI01", COHORT_TUMOR_FRACTION_COLUMN].iloc[0] == pytest.approx(0.2)
    assert "fragment_size_median" in matrix.columns
    assert "short_long_ratio_window_1" in matrix.columns
    assert "cnv_logr_bin_5" in matrix.columns


def test_arm_zscores_are_keyed_by_arm_not_by_per_sample_index(tmp_path: Path) -> None:
    """The regression this module exists to prevent.

    A sample's arm_zscores.csv holds only the arms it has windows on, so SNY02 (chr1
    only) writes two rows where SNY01 writes four. Pivoting on row position would put
    SNY02's *chr1* z-scores under the same columns as the others' -- which happens to be
    right here and would be wrong the moment the missing arm is not the last one. Keying
    on the `arm` label makes the absence explicit instead.
    """
    _build_cohort(tmp_path / "derived")
    assembler = CohortFeatureMatrixAssembler(_config(tmp_path / "derived", tmp_path / "cohort.csv"))

    matrix = assembler.assemble().set_index("sample_id")

    assert {"arm_zscore_chr1p", "arm_zscore_chr1q", "arm_zscore_chr2p", "arm_zscore_chr2q"} <= set(matrix.columns)
    assert matrix.loc["SNY02", ["arm_zscore_chr2p", "arm_zscore_chr2q"]].isna().all()
    assert matrix.loc["SNY01", ["arm_zscore_chr2p", "arm_zscore_chr2q"]].notna().all()

    # SNY02's own two arms land under the chr1 labels, matching its arm_zscores.csv.
    raw = read_csv(tmp_path / "derived" / "Snyder" / "SNY02" / "arm_zscores.csv").set_index("arm")["zscore"]
    assert matrix.loc["SNY02", "arm_zscore_chr1p"] == pytest.approx(raw.loc["chr1p"])
    assert matrix.loc["SNY02", "arm_zscore_chr1q"] == pytest.approx(raw.loc["chr1q"])


def test_missing_windows_become_nan_rather_than_shifting_columns(tmp_path: Path) -> None:
    _build_cohort(tmp_path / "derived")
    assembler = CohortFeatureMatrixAssembler(_config(tmp_path / "derived", tmp_path / "cohort.csv"))

    matrix = assembler.assemble().set_index("sample_id")

    # window_index 41 is chr2's first window: present for SNY01/CRI01, absent for SNY02.
    assert matrix.loc["SNY02", "short_long_ratio_window_41"] != matrix.loc["SNY02", "short_long_ratio_window_41"]
    assert matrix.loc["SNY01", "short_long_ratio_window_41"] == matrix.loc["SNY01", "short_long_ratio_window_41"]


def test_samples_without_a_completion_flag_are_skipped(tmp_path: Path) -> None:
    derived_root = tmp_path / "derived"
    _build_cohort(derived_root)
    _write_sample(derived_root, "Snyder", "SNY03", CHROM_WINDOW_COUNTS, tumor_fraction=0.3, random_state=3, complete=False)

    assembler = CohortFeatureMatrixAssembler(_config(derived_root, tmp_path / "cohort.csv"))

    assert "SNY03" not in set(assembler.assemble()["sample_id"])


def test_duplicate_sample_id_across_datasets_is_rejected(tmp_path: Path) -> None:
    derived_root = tmp_path / "derived"
    _build_cohort(derived_root)
    _write_sample(derived_root, "Cristiano", "SNY01", CHROM_WINDOW_COUNTS, tumor_fraction=0.4, random_state=4)

    assembler = CohortFeatureMatrixAssembler(_config(derived_root, tmp_path / "cohort.csv"))

    with pytest.raises(ValueError, match="SNY01"):
        assembler.assemble()


def test_write_and_summarise(tmp_path: Path) -> None:
    _build_cohort(tmp_path / "derived")
    assembler = CohortFeatureMatrixAssembler(_config(tmp_path / "derived", tmp_path / "out" / "cohort.csv"))

    matrix = assembler.assemble()
    output_path = assembler.write(matrix)
    summary = assembler.summarise(matrix)

    assert read_csv(output_path).shape == matrix.shape
    assert summary["n_samples"] == 3
    assert summary["samples_by_source_dataset"] == {"Snyder": 2, "Cristiano": 1}
    # SNY02 lacks chr2's windows and arms, so those columns are not complete cohort-wide.
    assert summary["complete_feature_columns"] < summary["n_feature_columns"]
    assert summary["tumor_fraction_nonzero"] == 2
