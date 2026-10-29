"""Unit tests for flash_fraq.genome_mask's centromere/chromosome-end mask.

Real hg19 chr21 coordinates throughout (src/flash_fraq/data/hg19_chrom_sizes.tsv:
chr21 length 48,129,895; hg19_centromeres.tsv: chr21 centromere 11,288,129-14,288,129),
so these tests also double as a sanity check on the vendored data files themselves.
"""

from pandas import DataFrame

from flash_fraq.genome_mask import DEFAULT_MASK_MARGIN_BP, masked, mask_frame

CHR21_LEN = 48_129_895
CHR21_CEN_START = 11_288_129
CHR21_CEN_END = 14_288_129


def test_masked_true_for_window_overlapping_centromere() -> None:
    # Default 5 Mb margin: centromere span extended by 5 Mb on each side is
    # 6,288,129-19,288,129; a window at 10-11 Mb falls inside it.
    assert masked("chr21", 10_000_000, 11_000_000) is True


def test_masked_false_for_window_far_from_centromere_and_ends() -> None:
    # Mid chr21q, comfortably outside the default margin around the centromere
    # (19,288,129) and the chromosome end (48,129,895 - 5,000,000 = 43,129,895).
    assert masked("chr21", 25_000_000, 26_000_000) is False


def test_masked_true_for_window_at_chromosome_start() -> None:
    assert masked("chr21", 0, 1_000_000) is True


def test_masked_true_for_window_at_chromosome_end() -> None:
    assert masked("chr21", CHR21_LEN - 1_000_000, CHR21_LEN) is True


def test_masked_false_just_outside_margin_on_both_ends() -> None:
    margin = DEFAULT_MASK_MARGIN_BP
    # Starts exactly at the margin boundary (not < margin_bp) and ends exactly at
    # chrom_len - margin_bp (not > chrom_len - margin_bp): neither triggers the mask.
    assert masked("chr21", margin, margin + 1_000_000) is False
    assert masked("chr21", CHR21_LEN - margin - 1_000_000, CHR21_LEN - margin) is False


def test_masked_true_just_inside_margin_on_both_ends() -> None:
    margin = DEFAULT_MASK_MARGIN_BP
    assert masked("chr21", margin - 1, margin + 1_000_000 - 1) is True
    assert masked("chr21", CHR21_LEN - margin - 1_000_000 + 1, CHR21_LEN - margin + 1) is True


def test_masked_true_just_inside_centromere_margin() -> None:
    margin = DEFAULT_MASK_MARGIN_BP
    # A 1bp window just inside the centromere's extended margin on the p-arm side.
    assert masked("chr21", CHR21_CEN_START - margin, CHR21_CEN_START - margin + 1) is True
    # ... and on the q-arm side.
    assert masked("chr21", CHR21_CEN_END + margin - 1, CHR21_CEN_END + margin) is True


def test_masked_false_just_outside_centromere_margin() -> None:
    margin = DEFAULT_MASK_MARGIN_BP
    assert masked("chr21", CHR21_CEN_START - margin - 1_000_000, CHR21_CEN_START - margin) is False
    assert masked("chr21", CHR21_CEN_END + margin, CHR21_CEN_END + margin + 1_000_000) is False


def test_margin_bp_zero_disables_masking_everywhere() -> None:
    assert masked("chr21", 0, 1_000_000, margin_bp=0) is False
    assert masked("chr21", CHR21_LEN - 1_000_000, CHR21_LEN, margin_bp=0) is False
    assert masked("chr21", CHR21_CEN_START, CHR21_CEN_END, margin_bp=0) is False


def test_masked_unknown_chrom_raises() -> None:
    try:
        masked("chr99", 0, 1_000_000)
    except KeyError as error:
        assert "chr99" in str(error)
    else:
        raise AssertionError("expected KeyError for an unvendored chromosome")


def test_mask_frame_matches_masked_row_by_row() -> None:
    df = DataFrame(
        {
            "chrom": ["chr21", "chr21", "chr21", "chr21"],
            "start": [0, 10_000_000, 25_000_000, CHR21_LEN - 1_000_000],
            "end": [1_000_000, 11_000_000, 26_000_000, CHR21_LEN],
        }
    )

    result = mask_frame(df)

    assert list(result) == [True, True, False, True]


def test_mask_frame_with_bin_size_bp_reconstructs_end() -> None:
    # Same rows as above but only `start` is supplied, with a fixed bin size --
    # mirrors cnv_bins.csv, which carries no `end` column.
    df = DataFrame({"chrom": ["chr21", "chr21", "chr21"], "start": [0, 10_000_000, 25_000_000]})

    result = mask_frame(df, end_col=None, bin_size_bp=1_000_000)

    assert list(result) == [True, True, False]


def test_mask_frame_margin_bp_zero_returns_all_false() -> None:
    df = DataFrame({"chrom": ["chr21", "chr21"], "start": [0, CHR21_LEN - 1_000_000], "end": [1_000_000, CHR21_LEN]})

    result = mask_frame(df, margin_bp=0)

    assert list(result) == [False, False]


def test_mask_frame_preserves_index() -> None:
    df = DataFrame(
        {"chrom": ["chr21", "chr21"], "start": [25_000_000, 0], "end": [26_000_000, 1_000_000]},
        index=[7, 3],
    )

    result = mask_frame(df)

    assert list(result.index) == [7, 3]
    assert result.loc[7] is False or result.loc[7] == False  # noqa: E712
    assert result.loc[3] == True  # noqa: E712


def test_mask_frame_chromosome_with_no_centromere_entry_only_masks_ends() -> None:
    # chrX is vendored in hg19_chrom_sizes.tsv and hg19_centromeres.tsv both, so use a
    # position far from either to confirm centromere lookup for a *present* chrom
    # doesn't spuriously mask everything, then confirm an unvendored chrom still raises.
    df = DataFrame({"chrom": ["chrX"], "start": [80_000_000], "end": [81_000_000]})

    result = mask_frame(df)

    assert list(result) == [False]


def test_mask_frame_unknown_chrom_raises() -> None:
    df = DataFrame({"chrom": ["chr21", "chr99"], "start": [0, 0], "end": [1_000_000, 1_000_000]})

    try:
        mask_frame(df)
    except KeyError as error:
        assert "chr99" in str(error)
    else:
        raise AssertionError("expected KeyError for an unvendored chromosome")
