from pathlib import Path

import numpy
import pytest
from pandas import DataFrame, read_csv

from flash_fraq.gc_correction import HG19_CENTROMERE_MB, GcCorrectionConfig, GcCorrector, arm_id


def _make_synthetic_windows(n_windows: int = 200) -> DataFrame:
    # Shuffled (not genomic-order-monotonic) so the GC-correlated trend and the
    # arm-level step below are independent signals, as they would be in real data.
    gc_content = numpy.random.default_rng(0).permutation(numpy.linspace(0.30, 0.65, n_windows))
    true_short = 1000.0
    true_long = 2000.0
    short_bias = 900.0 * (gc_content - gc_content.mean())
    long_bias = -500.0 * (gc_content - gc_content.mean())
    # Windows spread across chr1/chr2 so arm bucketing has more than one arm to work with.
    chrom = ["chr1"] * (n_windows // 2) + ["chr2"] * (n_windows - n_windows // 2)
    start = [i * 5_000_000 for i in range(n_windows // 2)] + [i * 5_000_000 for i in range(n_windows - n_windows // 2)]
    # A GC-independent per-arm offset so arm-level means differ meaningfully after GC
    # correction flattens the GC-correlated trend (otherwise all arms converge to the
    # same corrected ratio and z-scoring divides by a near-zero standard deviation).
    is_p_arm = [s < HG19_CENTROMERE_MB[c] * 1_000_000 for c, s in zip(chrom, start)]
    arm_offset = numpy.array([40.0 if p else -30.0 for p in is_p_arm])
    return DataFrame(
        {
            "sample_id": "test_sample",
            "window_index": range(1, n_windows + 1),
            "chrom": chrom,
            "start": start,
            "gc_content": gc_content,
            "short_count": true_short + short_bias + arm_offset,
            "long_count": true_long + long_bias,
        }
    )


def test_correct_removes_gc_correlated_bias(tmp_path: Path) -> None:
    windows = _make_synthetic_windows()
    input_path = tmp_path / "windows_raw.tsv"
    output_path = tmp_path / "windows_corrected.tsv"
    windows.to_csv(input_path, sep="\t", index=False)

    config = GcCorrectionConfig.from_sources(
        input_windows_path=input_path,
        output_windows_path=output_path,
        span=0.75,
        mask_margin_bp=0,
    )
    corrector = GcCorrector(config)
    corrected = corrector.run()

    raw_ratio = windows["short_count"] / windows["long_count"]
    raw_correlation = raw_ratio.corr(windows["gc_content"])
    corrected_correlation = corrected["short_long_ratio"].corr(corrected["gc_content"])

    assert abs(raw_correlation) > 0.5
    assert abs(corrected_correlation) < 0.1
    assert output_path.exists()


def test_correct_counts_preserve_genome_wide_median(tmp_path: Path) -> None:
    windows = _make_synthetic_windows()
    config = GcCorrectionConfig.from_sources(
        input_windows_path=tmp_path / "unused_raw.tsv",
        output_windows_path=tmp_path / "unused_corrected.tsv",
        mask_margin_bp=0,
    )
    corrector = GcCorrector(config)
    corrected = corrector.correct(windows)

    assert corrected["short_count_gc_corrected"].median() == pytest.approx(
        windows["short_count"].median(), abs=1.0
    )


def test_write_window_ratios_matches_feature_matrix_schema(tmp_path: Path) -> None:
    windows = _make_synthetic_windows()
    config = GcCorrectionConfig.from_sources(
        input_windows_path=tmp_path / "windows_raw.tsv",
        output_windows_path=tmp_path / "windows_corrected.tsv",
        mask_margin_bp=0,
    )
    windows.to_csv(config.input_windows_path, sep="\t", index=False)
    corrector = GcCorrector(config)
    corrected = corrector.run()

    window_ratios_path = corrector.write_window_ratios(corrected)

    written = read_csv(window_ratios_path)
    assert list(written.columns) == [
        "sample_id", "window_index", "short_long_ratio", "gc_corrected_total_coverage", "gc_corrected_short_coverage",
    ]
    assert written.shape[0] == windows.shape[0]


def test_compute_arm_zscores_buckets_windows_by_chromosome_arm(tmp_path: Path) -> None:
    windows = _make_synthetic_windows()
    config = GcCorrectionConfig.from_sources(
        input_windows_path=tmp_path / "windows_raw.tsv",
        output_windows_path=tmp_path / "windows_corrected.tsv",
        mask_margin_bp=0,
    )
    corrector = GcCorrector(config)
    corrected = corrector.correct(windows)

    arm_zscores = corrector.compute_arm_zscores(corrected)

    assert set(arm_zscores["sample_id"]) == {"test_sample"}
    assert arm_zscores["zscore"].mean() == pytest.approx(0.0, abs=1e-6)
    assert arm_zscores.shape[0] >= 2
    assert arm_zscores["zscore"].std() > 0.1
    assert list(arm_zscores["arm"]) == sorted(
        {arm_id(chrom, start) for chrom, start in zip(windows["chrom"], windows["start"])}
    )


def test_write_arm_zscores_carries_the_arm_label(tmp_path: Path) -> None:
    """The label, not the row's position, is the cross-sample key.

    A sample only has windows on the arms its coverage reaches, so the arm set is ragged
    across a cohort and a positional index would put different arms in one column of any
    join -- the defect that cost a 547-sample back-fill.
    """
    windows = _make_synthetic_windows()
    config = GcCorrectionConfig.from_sources(
        input_windows_path=tmp_path / "windows_raw.tsv",
        output_windows_path=tmp_path / "windows_corrected.tsv",
        mask_margin_bp=0,
    )
    windows.to_csv(config.input_windows_path, sep="\t", index=False)
    corrector = GcCorrector(config)
    corrected = corrector.run()

    arm_zscore_path = corrector.write_arm_zscores(corrected)

    written = read_csv(arm_zscore_path)
    assert list(written.columns) == ["sample_id", "arm", "zscore"]
    assert list(written["arm"]) == sorted(
        {arm_id(chrom, start) for chrom, start in zip(windows["chrom"], windows["start"])}
    )
    assert written["zscore"].to_numpy() == pytest.approx(
        corrector.compute_arm_zscores(corrected)["zscore"].to_numpy()
    )


def _make_chr21_windows() -> DataFrame:
    """Real hg19 chr21 5 Mb windows (chr21 is acrocentric: hg19.fa.fai length
    48,129,895; UCSC gap-table centromere 11,288,129-14,288,129), covering the whole
    chromosome so a masking test can check both the chromosome-start/end margin and
    the centromere margin, and the acrocentric-arm-fully-masked case, against real
    coordinates rather than made-up ones.
    """
    chrom_len = 48_129_895
    starts = list(range(0, chrom_len, 5_000_000))
    ends = [min(start + 5_000_000, chrom_len) for start in starts]
    n = len(starts)
    gc_content = numpy.linspace(0.30, 0.65, n)
    return DataFrame(
        {
            "sample_id": "test_sample",
            "window_index": range(1, n + 1),
            "chrom": "chr21",
            "start": starts,
            "end": ends,
            "gc_content": gc_content,
            "short_count": 1000.0 + 100.0 * gc_content,
            "long_count": 2000.0 - 50.0 * gc_content,
        }
    )


# Manually verified against src/flash_fraq/data/hg19_centromeres.tsv's chr21 span
# (11,288,129-14,288,129) and hg19_chrom_sizes.tsv's chr21 length (48,129,895) at the
# default 5 Mb margin: starts 0/5e6/10e6/15e6 fall within margin_bp of the centromere
# (chr21p is entirely within it -- the acrocentric case), 40e6/45e6 fall within
# margin_bp of the chromosome end, and only 20e6/25e6/30e6/35e6 survive, all on chr21q.
CHR21_MASKED_STARTS = {0, 5_000_000, 10_000_000, 15_000_000, 40_000_000, 45_000_000}
CHR21_UNMASKED_STARTS = {20_000_000, 25_000_000, 30_000_000, 35_000_000}


def test_correct_drops_masked_windows_before_the_loess_fit(tmp_path: Path) -> None:
    windows = _make_chr21_windows()
    config = GcCorrectionConfig.from_sources(
        input_windows_path=tmp_path / "windows_raw.tsv",
        output_windows_path=tmp_path / "windows_corrected.tsv",
    )
    corrector = GcCorrector(config)

    corrected = corrector.correct(windows)

    assert set(corrected["start"]) == CHR21_UNMASKED_STARTS
    assert corrected.shape[0] == len(CHR21_UNMASKED_STARTS)


def test_mask_margin_bp_zero_disables_masking(tmp_path: Path) -> None:
    windows = _make_chr21_windows()
    config = GcCorrectionConfig.from_sources(
        input_windows_path=tmp_path / "windows_raw.tsv",
        output_windows_path=tmp_path / "windows_corrected.tsv",
        mask_margin_bp=0,
    )
    corrector = GcCorrector(config)

    corrected = corrector.correct(windows)

    assert set(corrected["start"]) == set(windows["start"])


def test_write_window_ratios_excludes_masked_windows(tmp_path: Path) -> None:
    windows = _make_chr21_windows()
    config = GcCorrectionConfig.from_sources(
        input_windows_path=tmp_path / "windows_raw.tsv",
        output_windows_path=tmp_path / "windows_corrected.tsv",
    )
    windows.to_csv(config.input_windows_path, sep="\t", index=False)
    corrector = GcCorrector(config)
    corrected = corrector.run()

    window_ratios_path = corrector.write_window_ratios(corrected)
    written = read_csv(window_ratios_path)

    masked_window_indices = set(
        windows.loc[windows["start"].isin(CHR21_MASKED_STARTS), "window_index"]
    )
    assert not masked_window_indices & set(written["window_index"])
    assert written.shape[0] == len(CHR21_UNMASKED_STARTS)


def test_acrocentric_p_arm_fully_masked_is_dropped_not_nan(tmp_path: Path) -> None:
    """chr21p sits entirely inside the centromere margin (the acrocentric case): every
    one of its windows is masked, so it must not appear in arm_zscores.csv at all --
    consistent with how a sample's naturally ragged per-arm coverage already makes
    compute_arm_zscores/pivot_arm_zscores treat a missing arm (dropped, never a NaN row
    or column)."""
    windows = _make_chr21_windows()
    config = GcCorrectionConfig.from_sources(
        input_windows_path=tmp_path / "windows_raw.tsv",
        output_windows_path=tmp_path / "windows_corrected.tsv",
    )
    corrector = GcCorrector(config)
    corrected = corrector.correct(windows)

    arm_zscores = corrector.compute_arm_zscores(corrected)

    assert list(arm_zscores["arm"]) == ["chr21q"]
    assert "chr21p" not in set(arm_zscores["arm"])
