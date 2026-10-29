from pathlib import Path

import pysam
import pytest

from flash_fraq.tumor_batch_sampling import TumorBatchSampler, TumorDownsampleTarget, _GOLDEN_RATIO_CONJUGATE

CONTIG = "chr1"
CONTIG_LENGTH = 1_000_000
N_READS = 300


def _build_source_bam(tmp_path: Path, n_reads: int = N_READS) -> Path:
    bam_path = tmp_path / "tumor_source.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i in range(n_reads):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = f"read_{i}"
            read.query_sequence = "A" * 50
            read.flag = 0
            read.reference_id = 0
            read.reference_start = i
            read.mapping_quality = 60
            read.cigar = [(0, 50)]
            bam.write(read)
    return bam_path


def _read_names(bam_path: Path) -> list[str]:
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        return [read.query_name for read in bam]


@pytest.fixture
def source_bam(tmp_path: Path) -> Path:
    return _build_source_bam(tmp_path)


@pytest.fixture
def reference_fasta(tmp_path: Path) -> Path:
    # BAM (not CRAM) decoding does not touch this file's contents; it only
    # needs to exist for run()'s reference_filename argument.
    path = tmp_path / "ref.fa"
    path.write_text(f">{CONTIG}\n{'A' * CONTIG_LENGTH}\n")
    return path


def test_target_rejects_fraction_outside_zero_one(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="fraction must be in"):
        TumorDownsampleTarget(output_bam=tmp_path / "out.bam", seed=1, fraction=0.0)
    with pytest.raises(ValueError, match="fraction must be in"):
        TumorDownsampleTarget(output_bam=tmp_path / "out.bam", seed=1, fraction=1.5)


def test_target_phase_matches_the_weyl_sequence_formula(tmp_path: Path) -> None:
    target = TumorDownsampleTarget(output_bam=tmp_path / "out.bam", seed=7, fraction=0.5)
    assert target.phase == pytest.approx((7 * _GOLDEN_RATIO_CONJUGATE) % 1.0)


def test_run_with_no_targets_returns_zero_and_writes_nothing(source_bam: Path, reference_fasta: Path) -> None:
    sampler = TumorBatchSampler()
    total_reads = sampler.run(source_bam, reference_fasta, [])
    assert total_reads == 0


def test_run_returns_the_total_source_read_count(
    source_bam: Path, reference_fasta: Path, tmp_path: Path
) -> None:
    sampler = TumorBatchSampler()
    target = TumorDownsampleTarget(output_bam=tmp_path / "out" / "row.tumor_sub.bam", seed=1, fraction=0.5)
    total_reads = sampler.run(source_bam, reference_fasta, [target])
    assert total_reads == N_READS


def test_run_with_fraction_one_retains_every_read(
    source_bam: Path, reference_fasta: Path, tmp_path: Path
) -> None:
    sampler = TumorBatchSampler()
    target = TumorDownsampleTarget(output_bam=tmp_path / "out" / "row.tumor_sub.bam", seed=3, fraction=1.0)
    sampler.run(source_bam, reference_fasta, [target])

    assert set(_read_names(target.output_bam)) == set(_read_names(source_bam))


def test_run_is_reproducible_for_the_same_seed_and_fraction(
    source_bam: Path, reference_fasta: Path, tmp_path: Path
) -> None:
    sampler = TumorBatchSampler()
    target_a = TumorDownsampleTarget(output_bam=tmp_path / "a" / "row.tumor_sub.bam", seed=11, fraction=0.3)
    target_b = TumorDownsampleTarget(output_bam=tmp_path / "b" / "row.tumor_sub.bam", seed=11, fraction=0.3)

    sampler.run(source_bam, reference_fasta, [target_a])
    sampler.run(source_bam, reference_fasta, [target_b])

    names_a = _read_names(target_a.output_bam)
    names_b = _read_names(target_b.output_bam)
    assert names_a  # sanity: a non-trivial fraction actually retained something
    assert names_a == names_b


def test_run_produces_distinct_subsets_for_distinct_seeds_at_the_same_fraction(
    source_bam: Path, reference_fasta: Path, tmp_path: Path
) -> None:
    sampler = TumorBatchSampler()
    target_a = TumorDownsampleTarget(output_bam=tmp_path / "a" / "row.tumor_sub.bam", seed=1, fraction=0.5)
    target_b = TumorDownsampleTarget(output_bam=tmp_path / "b" / "row.tumor_sub.bam", seed=2, fraction=0.5)

    sampler.run(source_bam, reference_fasta, [target_a, target_b])

    names_a = set(_read_names(target_a.output_bam))
    names_b = set(_read_names(target_b.output_bam))
    assert names_a != names_b


def test_run_writes_all_targets_from_a_single_pass_matching_independent_single_target_runs(
    source_bam: Path, reference_fasta: Path, tmp_path: Path
) -> None:
    """The whole point of the batch sampler: running N targets together in
    one call must give exactly the same per-target result as running each
    target alone -- the batching (and the numpy-vectorized decision matrix)
    must not change any individual target's selection."""
    targets = [
        TumorDownsampleTarget(output_bam=tmp_path / "batched" / f"row{i}.tumor_sub.bam", seed=i, fraction=0.4)
        for i in range(5)
    ]
    TumorBatchSampler(batch_size=37).run(source_bam, reference_fasta, targets)

    for i, target in enumerate(targets):
        solo_target = TumorDownsampleTarget(
            output_bam=tmp_path / "solo" / f"row{i}.tumor_sub.bam", seed=i, fraction=0.4
        )
        TumorBatchSampler().run(source_bam, reference_fasta, [solo_target])
        assert set(_read_names(target.output_bam)) == set(_read_names(solo_target.output_bam))


def test_run_creates_missing_output_parent_directories(
    source_bam: Path, reference_fasta: Path, tmp_path: Path
) -> None:
    sampler = TumorBatchSampler()
    target = TumorDownsampleTarget(
        output_bam=tmp_path / "nested" / "does" / "not" / "exist" / "row.tumor_sub.bam", seed=1, fraction=0.5
    )
    sampler.run(source_bam, reference_fasta, [target])
    assert target.output_bam.exists()


def test_both_mates_of_a_pair_get_the_same_inclusion_decision(
    tmp_path: Path, reference_fasta: Path
) -> None:
    bam_path = tmp_path / "paired.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for i in range(100):
            for mate_flag in (0x40, 0x80):
                read = pysam.AlignedSegment(bam.header)
                read.query_name = f"pair_{i}"
                read.query_sequence = "A" * 50
                read.flag = 0x1 | mate_flag
                read.reference_id = 0
                read.reference_start = i
                read.mapping_quality = 60
                read.cigar = [(0, 50)]
                bam.write(read)

    sampler = TumorBatchSampler()
    target = TumorDownsampleTarget(output_bam=tmp_path / "out" / "row.tumor_sub.bam", seed=5, fraction=0.5)
    sampler.run(bam_path, reference_fasta, [target])

    names = _read_names(target.output_bam)
    # Every retained query_name should appear exactly twice (both mates kept
    # together) -- never a lone mate.
    from collections import Counter

    counts = Counter(names)
    assert counts
    assert set(counts.values()) == {2}
