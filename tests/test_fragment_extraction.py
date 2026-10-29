from pathlib import Path

import pysam
import pytest

from flash_fraq.fragment_extraction import (
    FRAGMENT_BED_INTERMEDIATES,
    FragmentExtractionConfig,
    FragmentExtractor,
    cleanup_fragment_intermediates,
)


PAIRED = 0x1
PROPER_PAIR = 0x2
READ1 = 0x40
READ2 = 0x80
MATE_REVERSE = 0x20
REVERSE = 0x10
DUPLICATE = 0x400
SECONDARY = 0x100
SUPPLEMENTARY = 0x800


def _add_pair(
    bam: pysam.AlignmentFile,
    name: str,
    read1_start: int,
    read2_start: int,
    fragment_len: int,
    extra_read1_flags: int = 0,
) -> None:
    read1 = pysam.AlignedSegment(bam.header)
    read1.query_name = name
    read1.query_sequence = "A" * 50
    read1.flag = PAIRED | PROPER_PAIR | READ1 | MATE_REVERSE | extra_read1_flags
    read1.reference_id = 0
    read1.reference_start = read1_start
    read1.mapping_quality = 60
    read1.cigar = [(0, 50)]
    read1.next_reference_id = 0
    read1.next_reference_start = read2_start
    read1.template_length = fragment_len
    bam.write(read1)

    read2 = pysam.AlignedSegment(bam.header)
    read2.query_name = name
    read2.query_sequence = "A" * 50
    read2.flag = PAIRED | PROPER_PAIR | READ2 | REVERSE
    read2.reference_id = 0
    read2.reference_start = read2_start
    read2.mapping_quality = 60
    read2.cigar = [(0, 50)]
    read2.next_reference_id = 0
    read2.next_reference_start = read1_start
    read2.template_length = -fragment_len
    bam.write(read2)


@pytest.fixture
def synthetic_bam(tmp_path: Path) -> Path:
    bam_path = tmp_path / "synthetic.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": "chr1", "LN": 10_000_000}],
    }
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        # Kept: fragment_start=100, tlen=200 -> midpoint=200, in window 1 / bin 1.
        _add_pair(bam, "pair1_kept", read1_start=100, read2_start=250, fragment_len=200)
        # Kept: fragment_start=6_000_100, tlen=150 -> midpoint=6_000_175, in window 2 / bin 7.
        _add_pair(bam, "pair2_kept", read1_start=6_000_100, read2_start=6_000_200, fragment_len=150)
        # Excluded: duplicate.
        _add_pair(
            bam, "pair3_duplicate", read1_start=300, read2_start=450, fragment_len=200,
            extra_read1_flags=DUPLICATE,
        )
        # Excluded: secondary alignment.
        _add_pair(
            bam, "pair4_secondary", read1_start=400, read2_start=550, fragment_len=200,
            extra_read1_flags=SECONDARY,
        )
        # Excluded: supplementary alignment.
        _add_pair(
            bam, "pair5_supplementary", read1_start=500, read2_start=650, fragment_len=200,
            extra_read1_flags=SUPPLEMENTARY,
        )
        # Excluded: not properly paired (flag omits PROPER_PAIR).
        read1 = pysam.AlignedSegment(bam.header)
        read1.query_name = "pair6_not_proper"
        read1.query_sequence = "A" * 50
        read1.flag = PAIRED | READ1 | MATE_REVERSE
        read1.reference_id = 0
        read1.reference_start = 700
        read1.mapping_quality = 60
        read1.cigar = [(0, 50)]
        read1.next_reference_id = 0
        read1.next_reference_start = 850
        read1.template_length = 200
        bam.write(read1)
    return bam_path


@pytest.fixture
def delfi_windows_bed(tmp_path: Path) -> Path:
    bed_path = tmp_path / "delfi_5mb_windows.bed"
    bed_path.write_text(
        "chr1\t0\t5000000\twindow1\n"
        "chr1\t5000000\t10000000\twindow2\n"
    )
    return bed_path


@pytest.fixture
def ichor_bins_bed(tmp_path: Path) -> Path:
    bed_path = tmp_path / "ichor_1mb_bins.bed"
    lines = [f"chr1\t{i * 1_000_000}\t{(i + 1) * 1_000_000}\tbin{i + 1}\n" for i in range(10)]
    bed_path.write_text("".join(lines))
    return bed_path


@pytest.fixture
def reference_fasta(tmp_path: Path) -> Path:
    import pysam as _pysam

    fasta_path = tmp_path / "ref.fa"
    # 50% GC throughout, so gc_content for both windows is deterministic and simple.
    sequence = "GCAT" * 2_500_000
    with fasta_path.open("w") as handle:
        handle.write(">chr1\n")
        for i in range(0, len(sequence), 70):
            handle.write(sequence[i:i + 70] + "\n")
    _pysam.faidx(str(fasta_path))
    return fasta_path


@pytest.fixture
def config(
    synthetic_bam: Path, delfi_windows_bed: Path, ichor_bins_bed: Path, reference_fasta: Path, tmp_path: Path,
) -> FragmentExtractionConfig:
    return FragmentExtractionConfig.from_sources(
        input_bam=synthetic_bam,
        output_dir=tmp_path / "output",
        delfi_windows_bed=delfi_windows_bed,
        ichor_bins_bed=ichor_bins_bed,
        reference_fasta=reference_fasta,
        sample_id="test_sample",
    )


def test_extract_midpoints_filters_and_computes_expected_values(config: FragmentExtractionConfig) -> None:
    extractor = FragmentExtractor(config)

    midpoints_bed = extractor.extract_midpoints()

    lines = midpoints_bed.read_text().strip().split("\n")
    assert lines == [
        "chr1\t200\t201\t200",
        "chr1\t6000175\t6000176\t150",
    ]


def test_intersect_windows_assigns_fragments_to_delfi_and_ichor_windows(config: FragmentExtractionConfig) -> None:
    extractor = FragmentExtractor(config)

    delfi_out, ichor_out = extractor.run()

    delfi_lines = delfi_out.read_text().strip().split("\n")
    assert delfi_lines == [
        "chr1\t200\t201\t200\tchr1\t0\t5000000\twindow1",
        "chr1\t6000175\t6000176\t150\tchr1\t5000000\t10000000\twindow2",
    ]

    ichor_lines = ichor_out.read_text().strip().split("\n")
    assert ichor_lines == [
        "chr1\t200\t201\t200\tchr1\t0\t1000000\tbin1",
        "chr1\t6000175\t6000176\t150\tchr1\t6000000\t7000000\tbin7",
    ]


def test_aggregate_window_counts_computes_short_long_counts_and_gc_content(
    config: FragmentExtractionConfig,
) -> None:
    extractor = FragmentExtractor(config)
    midpoints_bed = extractor.extract_midpoints()
    delfi_out, _ichor_out = extractor.intersect_windows(midpoints_bed)

    window_counts = extractor.aggregate_window_counts(delfi_out)

    assert list(window_counts["window_index"]) == [1, 2]
    assert list(window_counts["sample_id"]) == ["test_sample", "test_sample"]
    # fragment 1: tlen=200 -> long; fragment 2: tlen=150 -> short.
    assert list(window_counts["short_count"]) == [0, 1]
    assert list(window_counts["long_count"]) == [1, 0]
    assert window_counts["gc_content"].iloc[0] == pytest.approx(0.5)


def test_compute_fragment_stats_summarizes_tlen_distribution(config: FragmentExtractionConfig) -> None:
    extractor = FragmentExtractor(config)
    midpoints_bed = extractor.extract_midpoints()

    fragment_size_median, fragment_size_mode, fragment_size_entropy = extractor.compute_fragment_stats(midpoints_bed)

    assert fragment_size_median == pytest.approx(175)
    assert fragment_size_mode in (150, 200)
    assert fragment_size_entropy > 0


def test_run_full_writes_all_outputs(config: FragmentExtractionConfig) -> None:
    extractor = FragmentExtractor(config)

    delfi_out, ichor_out, window_counts_path, fragment_stats_path = extractor.run_full()

    assert delfi_out.exists()
    assert ichor_out.exists()
    assert window_counts_path.exists()
    assert fragment_stats_path.exists()

    fragment_stats_lines = fragment_stats_path.read_text().strip().split("\n")
    assert fragment_stats_lines[0] == "sample_id,fragment_size_median,fragment_size_mode,fragment_size_entropy"


def test_cleanup_fragment_intermediates_removes_the_three_beds(config: FragmentExtractionConfig) -> None:
    extractor = FragmentExtractor(config)
    delfi_out, ichor_out, window_counts_path, fragment_stats_path = extractor.run_full()
    assert delfi_out.exists() and ichor_out.exists()

    removed = cleanup_fragment_intermediates(config.output_dir)

    assert sorted(path.name for path in removed) == sorted(FRAGMENT_BED_INTERMEDIATES)
    assert not delfi_out.exists()
    assert not ichor_out.exists()
    assert not (config.output_dir / "fragment_midpoints.bed").exists()
    # Outputs that still have readers downstream are untouched.
    assert window_counts_path.exists()
    assert fragment_stats_path.exists()


def test_cleanup_fragment_intermediates_dry_run_does_not_delete(config: FragmentExtractionConfig) -> None:
    extractor = FragmentExtractor(config)
    delfi_out, ichor_out, _window_counts_path, _fragment_stats_path = extractor.run_full()

    removed = cleanup_fragment_intermediates(config.output_dir, dry_run=True)

    assert sorted(path.name for path in removed) == sorted(FRAGMENT_BED_INTERMEDIATES)
    assert delfi_out.exists()
    assert ichor_out.exists()


def test_cleanup_fragment_intermediates_is_idempotent_when_already_removed(tmp_path: Path) -> None:
    # No fragment-bed files ever written into this directory -- mirrors a retried
    # cleanup pass, or a driver run after an earlier pass already cleaned up.
    assert cleanup_fragment_intermediates(tmp_path) == []
