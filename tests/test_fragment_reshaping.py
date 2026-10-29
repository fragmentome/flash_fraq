"""Synthetic-BAM tests for in silico fragment-size reshaping
(docs/todo/completed/fragment-size-reshaping.md).

Two distinct fragment-size distributions stand in for the ticket's two real
populations: `_hcc1395_like_tlens` is a broad near-uniform spread (a sheared
cell-line library, imperfectly size-selected so it still carries some density
in the short range -- real mechanical shearing is not a hard cutoff) and
`_cfdna_like_tlens` is a tight Gaussian around 166bp (the canonical
mono-nucleosome cfDNA peak, Cristiano et al. 2019). The tests assert
`FragmentReshaper` moves the first population's median/mode toward the
second's, that retention is deterministic under a fixed seed, and that the
verification step's KS/distance check both fails on a deliberately
mismatched pair and passes once reshaping is applied -- the ticket's own
fallback for the "not trivially separable" negative control when a full
admixture-level classifier is out of scope.
"""

import random
from collections import Counter
from pathlib import Path

import pysam
import pytest

from flash_fraq.fragment_reshaping import (
    DEFAULT_MIN_MAPQ,
    FragmentReshaper,
    ReshapingConfig,
    acceptance_probabilities,
    compute_insert_size_histogram,
    ks_statistic_from_histograms,
    read_histogram_csv,
    reshape_bam,
    verify_reshaping,
    write_histogram_csv,
)

CONTIG = "chr1"
CONTIG_LENGTH = 400_000
READ_LENGTH = 50

PAIRED = 0x1
PROPER_PAIR = 0x2
READ1 = 0x40
READ2 = 0x80
MATE_REVERSE = 0x20
REVERSE = 0x10
DUPLICATE = 0x400


def _add_pair(
    bam: pysam.AlignmentFile,
    name: str,
    start: int,
    tlen: int,
    extra_flags: int = 0,
    mapq1: int = 60,
    mapq2: int = 60,
    set_mq_tag: bool = True,
    reference_id: int = 0,
) -> None:
    """`mapq1`/`mapq2` are each mate's own MAPQ; `set_mq_tag` (default True,
    matching real `bwa mem` output) writes each mate's `MQ` tag to the
    *other* mate's actual MAPQ, so `_meets_mapq_threshold`'s `min(own, MQ)`
    rule can be exercised against a realistic record. `reference_id`
    (default 0, i.e. the header's first `@SQ`) lets multi-contig fixtures
    place a pair on a contig other than the first."""
    mate_start = start + tlen - READ_LENGTH

    read1 = pysam.AlignedSegment(bam.header)
    read1.query_name = name
    read1.query_sequence = "A" * READ_LENGTH
    read1.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
    read1.flag = PAIRED | PROPER_PAIR | READ1 | MATE_REVERSE | extra_flags
    read1.reference_id = reference_id
    read1.reference_start = start
    read1.mapping_quality = mapq1
    read1.cigar = [(0, READ_LENGTH)]
    read1.next_reference_id = reference_id
    read1.next_reference_start = mate_start
    read1.template_length = tlen
    if set_mq_tag:
        read1.set_tag("MQ", mapq2)
    bam.write(read1)

    read2 = pysam.AlignedSegment(bam.header)
    read2.query_name = name
    read2.query_sequence = "A" * READ_LENGTH
    read2.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
    read2.flag = PAIRED | PROPER_PAIR | READ2 | REVERSE | extra_flags
    read2.reference_id = reference_id
    read2.reference_start = mate_start
    read2.mapping_quality = mapq2
    read2.cigar = [(0, READ_LENGTH)]
    read2.next_reference_id = reference_id
    read2.next_reference_start = start
    read2.template_length = -tlen
    if set_mq_tag:
        read2.set_tag("MQ", mapq1)
    bam.write(read2)


def _build_bam(tmp_path: Path, name: str, tlens: list[int]) -> Path:
    unsorted_path = tmp_path / f"{name}.unsorted.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]}
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=header) as bam:
        position = 1000
        for i, tlen in enumerate(tlens):
            _add_pair(bam, f"{name}_{i}", position, tlen)
            position += tlen + 20
    sorted_path = tmp_path / f"{name}.bam"
    pysam.sort("-o", str(sorted_path), str(unsorted_path))
    pysam.index(str(sorted_path))
    return sorted_path


def _hcc1395_like_tlens(n: int, seed: int = 1) -> list[int]:
    """Broad, roughly-uniform spread over 100-500bp -- a sheared cell-line
    library. Includes some density down in the short range, as a real
    Covaris-sheared prep's size-selection tail would."""
    rng = random.Random(seed)
    return [rng.randint(100, 500) for _ in range(n)]


def _cfdna_like_tlens(n: int, seed: int = 2) -> list[int]:
    """Tight Gaussian around the mono-nucleosome cfDNA peak (~166bp)."""
    rng = random.Random(seed)
    return [max(80, min(250, round(rng.gauss(166, 15)))) for _ in range(n)]


def _reference_fasta(tmp_path: Path) -> Path:
    fasta_path = tmp_path / "ref.fa"
    sequence = "ACGT" * (CONTIG_LENGTH // 4 + 1)
    with fasta_path.open("w") as handle:
        handle.write(f">{CONTIG}\n")
        for i in range(0, CONTIG_LENGTH, 70):
            handle.write(sequence[i:i + 70] + "\n")
    pysam.faidx(str(fasta_path))
    return fasta_path


@pytest.fixture(scope="module")
def hcc1395_tlens() -> list[int]:
    return _hcc1395_like_tlens(4000)


@pytest.fixture(scope="module")
def cfdna_tlens() -> list[int]:
    return _cfdna_like_tlens(4000)


@pytest.fixture
def reference_fasta(tmp_path: Path) -> Path:
    return _reference_fasta(tmp_path)


@pytest.fixture
def source_bam(tmp_path: Path, hcc1395_tlens: list[int]) -> Path:
    return _build_bam(tmp_path, "hcc1395_like", hcc1395_tlens)


@pytest.fixture
def target_bam(tmp_path: Path, cfdna_tlens: list[int]) -> Path:
    return _build_bam(tmp_path, "cfdna_like", cfdna_tlens)


# --- histogram computation -----------------------------------------------


def test_compute_insert_size_histogram_counts_each_fragment_once(reference_fasta: Path, tmp_path: Path) -> None:
    bam_path = _build_bam(tmp_path, "fixed", [150, 150, 300])

    histogram = compute_insert_size_histogram(bam_path, reference_fasta)

    assert histogram == Counter({150: 2, 300: 1})
    assert sum(histogram.values()) == 3  # not 6 -- read2 is not double-counted


def test_compute_insert_size_histogram_respects_range_and_excludes_duplicates(
    reference_fasta: Path, tmp_path: Path
) -> None:
    unsorted_path = tmp_path / "with_dup.unsorted.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]}
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=header) as bam:
        _add_pair(bam, "kept", 1000, 200)
        _add_pair(bam, "duplicate", 2000, 200, extra_flags=DUPLICATE)
        _add_pair(bam, "too_short", 3000, 10)
        _add_pair(bam, "too_long", 4000, 900)
    sorted_path = tmp_path / "with_dup.bam"
    pysam.sort("-o", str(sorted_path), str(unsorted_path))
    pysam.index(str(sorted_path))

    histogram = compute_insert_size_histogram(sorted_path, reference_fasta, min_insert_size=30, max_insert_size=700)

    assert histogram == Counter({200: 1})


# --- MAPQ filter -------------------------------------------------------


def test_default_min_mapq_is_30() -> None:
    assert DEFAULT_MIN_MAPQ == 30


def test_compute_insert_size_histogram_excludes_low_mapq_pairs(reference_fasta: Path, tmp_path: Path) -> None:
    unsorted_path = tmp_path / "mapq.unsorted.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]}
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=header) as bam:
        _add_pair(bam, "high_mapq", 1000, 200, mapq1=60, mapq2=60)
        _add_pair(bam, "low_mapq", 2000, 200, mapq1=10, mapq2=10)
    sorted_path = tmp_path / "mapq.bam"
    pysam.sort("-o", str(sorted_path), str(unsorted_path))
    pysam.index(str(sorted_path))

    default_histogram = compute_insert_size_histogram(sorted_path, reference_fasta)
    assert default_histogram == Counter({200: 1})  # only the high-MAPQ pair

    unfiltered_histogram = compute_insert_size_histogram(sorted_path, reference_fasta, min_mapq=0)
    assert unfiltered_histogram == Counter({200: 2})  # both pairs, threshold disabled


def test_min_mapq_threshold_is_configurable(reference_fasta: Path, tmp_path: Path) -> None:
    unsorted_path = tmp_path / "mapq_cfg.unsorted.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]}
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=header) as bam:
        _add_pair(bam, "mapq_40", 1000, 200, mapq1=40, mapq2=40)
    sorted_path = tmp_path / "mapq_cfg.bam"
    pysam.sort("-o", str(sorted_path), str(unsorted_path))
    pysam.index(str(sorted_path))

    assert sum(compute_insert_size_histogram(sorted_path, reference_fasta, min_mapq=30).values()) == 1
    assert sum(compute_insert_size_histogram(sorted_path, reference_fasta, min_mapq=50).values()) == 0


def test_mapq_filter_considers_both_mates_via_mq_tag(reference_fasta: Path, tmp_path: Path) -> None:
    """One mate mapped confidently, the other poorly (as real bwa mem output
    for a read straddling a repeat can be): the pair must be excluded
    regardless of which mate's record `_meets_mapq_threshold` is evaluated
    against, since `MQ` carries the *other* mate's MAPQ into each record."""
    unsorted_path = tmp_path / "mixed_mapq.unsorted.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]}
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=header) as bam:
        _add_pair(bam, "mixed", 1000, 200, mapq1=60, mapq2=5)
    sorted_path = tmp_path / "mixed_mapq.bam"
    pysam.sort("-o", str(sorted_path), str(unsorted_path))
    pysam.index(str(sorted_path))

    histogram = compute_insert_size_histogram(sorted_path, reference_fasta, min_mapq=30)
    assert sum(histogram.values()) == 0  # the read1-keyed histogram still excludes the pair


def test_reshape_bam_drops_low_mapq_tumor_pairs(reference_fasta: Path, tmp_path: Path) -> None:
    unsorted_path = tmp_path / "reshape_mapq.unsorted.bam"
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}]}
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=header) as bam:
        _add_pair(bam, "high_mapq", 1000, 200, mapq1=60, mapq2=60)
        _add_pair(bam, "low_mapq", 2000, 200, mapq1=10, mapq2=10)
    sorted_path = tmp_path / "reshape_mapq.bam"
    pysam.sort("-o", str(sorted_path), str(unsorted_path))
    pysam.index(str(sorted_path))

    probabilities = {200: 1.0}
    output_bam = tmp_path / "out.bam"
    reshape_bam(sorted_path, reference_fasta, output_bam, probabilities, seed=1)

    with pysam.AlignmentFile(str(output_bam), "rb", reference_filename=str(reference_fasta)) as bam:
        names = {read.query_name for read in bam.fetch(until_eof=True)}
    assert names == {"high_mapq"}

    # min_mapq=0 disables the gate entirely, so the low-MAPQ pair survives too.
    output_bam_unfiltered = tmp_path / "out_unfiltered.bam"
    reshape_bam(sorted_path, reference_fasta, output_bam_unfiltered, probabilities, seed=1, min_mapq=0)
    with pysam.AlignmentFile(str(output_bam_unfiltered), "rb", reference_filename=str(reference_fasta)) as bam:
        names_unfiltered = {read.query_name for read in bam.fetch(until_eof=True)}
    assert names_unfiltered == {"high_mapq", "low_mapq"}


# --- acceptance probabilities ----------------------------------------------


def test_acceptance_probabilities_are_bounded_and_favor_target_overlap(
    reference_fasta: Path, source_bam: Path, target_bam: Path
) -> None:
    source_histogram = compute_insert_size_histogram(source_bam, reference_fasta)
    target_histogram = compute_insert_size_histogram(target_bam, reference_fasta)

    probabilities = acceptance_probabilities(source_histogram, target_histogram)

    assert probabilities  # non-empty
    assert all(0.0 <= p <= 1.0 for p in probabilities.values())
    # At least one insert size hits the ceiling -- rejection sampling always
    # scales so the best-covered bin is fully retained.
    assert max(probabilities.values()) == pytest.approx(1.0)
    # A size near the cfDNA peak (166bp) should be far more likely to survive
    # than a size out at the sheared-library tail (480bp), since the target
    # has essentially no mass there.
    near_peak = probabilities.get(166, probabilities.get(165, 0.0))
    far_tail = probabilities.get(480, 0.0)
    assert near_peak > far_tail


def test_acceptance_probabilities_rejects_non_overlapping_histograms() -> None:
    source_histogram = Counter({500: 10})
    target_histogram = Counter({100: 10})

    with pytest.raises(ValueError, match="no overlapping bins"):
        acceptance_probabilities(source_histogram, target_histogram)


def test_acceptance_probabilities_rejects_empty_histograms() -> None:
    with pytest.raises(ValueError, match="empty"):
        acceptance_probabilities(Counter(), Counter({100: 1}))
    with pytest.raises(ValueError, match="empty"):
        acceptance_probabilities(Counter({100: 1}), Counter())


# --- KS / histogram-distance verification ----------------------------------


def test_ks_statistic_is_zero_for_identical_histograms() -> None:
    histogram = Counter({100: 5, 150: 10, 200: 5})
    assert ks_statistic_from_histograms(histogram, histogram) == pytest.approx(0.0)


def test_ks_statistic_is_large_for_disjoint_histograms() -> None:
    statistic = ks_statistic_from_histograms(Counter({100: 10}), Counter({500: 10}))
    assert statistic == pytest.approx(1.0)


def test_verify_reshaping_fails_when_distributions_are_far_apart() -> None:
    reshaped = Counter({400: 50, 450: 50})
    target = Counter({150: 50, 160: 50})

    result = verify_reshaping(reshaped, target, threshold=0.15)

    assert result.passed is False
    assert result.statistic > result.threshold


# --- end-to-end reshaping ----------------------------------------------


def test_reshaper_moves_median_and_mode_toward_target_and_passes_verification(
    tmp_path: Path, reference_fasta: Path, source_bam: Path, target_bam: Path, cfdna_tlens: list[int]
) -> None:
    config = ReshapingConfig.from_sources(
        source_bam=source_bam,
        reference_fasta=reference_fasta,
        output_dir=tmp_path / "out",
        target_bam=target_bam,
        ks_threshold=0.15,
        seed=42,
        sample_name="hcc1395_reshaped",
    )
    reshaper = FragmentReshaper(config)

    source_histogram = reshaper.compute_source_histogram()
    result_summary = reshaper.run()

    verification = result_summary["verification"]
    assert result_summary["output_bam"].exists()
    assert result_summary["n_reshaped_fragments"] < result_summary["n_source_fragments"]

    # Moved toward the target: the reshaped median/mode land far closer to
    # the cfDNA-like ~166bp peak than the source's ~300bp uniform-range center.
    target_median = sorted(cfdna_tlens)[len(cfdna_tlens) // 2]
    source_median = compute_insert_size_histogram(source_bam, reference_fasta)
    from flash_fraq.fragment_reshaping import _histogram_median

    source_median_value = _histogram_median(source_median)
    assert abs(verification.reshaped_median - target_median) < abs(source_median_value - target_median)
    assert abs(verification.reshaped_median - target_median) < 30

    assert verification.passed is True
    assert verification.statistic <= config.ks_threshold
    assert result_summary["verification_report_path"].exists()

    # Histograms are written for later inspection.
    assert (config.output_dir / f"{config.sample_name}.source_histogram.csv").exists()
    assert (config.output_dir / f"{config.sample_name}.target_histogram.csv").exists()
    assert (config.output_dir / f"{config.sample_name}.reshaped_histogram.csv").exists()
    assert sum(source_histogram.values()) == result_summary["n_source_fragments"]


def test_reshape_is_deterministic_under_fixed_seed(
    tmp_path: Path, reference_fasta: Path, source_bam: Path, target_bam: Path
) -> None:
    source_histogram = compute_insert_size_histogram(source_bam, reference_fasta)
    target_histogram = compute_insert_size_histogram(target_bam, reference_fasta)
    probabilities = acceptance_probabilities(source_histogram, target_histogram)

    output_a = tmp_path / "a.bam"
    output_b = tmp_path / "b.bam"
    reshape_bam(source_bam, reference_fasta, output_a, probabilities, seed=7)
    reshape_bam(source_bam, reference_fasta, output_b, probabilities, seed=7)

    def _qnames(bam_path: Path) -> list[str]:
        with pysam.AlignmentFile(str(bam_path), "rb", reference_filename=str(reference_fasta)) as bam:
            return sorted(read.query_name for read in bam.fetch(until_eof=True))

    names_a = _qnames(output_a)
    names_b = _qnames(output_b)
    assert names_a == names_b
    assert len(names_a) > 0
    # Every retained qname appears exactly twice (both mates kept together).
    assert all(count == 2 for count in Counter(names_a).values())


def test_reshape_seed_changes_the_retained_set(reference_fasta: Path, source_bam: Path, target_bam: Path, tmp_path: Path) -> None:
    source_histogram = compute_insert_size_histogram(source_bam, reference_fasta)
    target_histogram = compute_insert_size_histogram(target_bam, reference_fasta)
    probabilities = acceptance_probabilities(source_histogram, target_histogram)

    output_a = tmp_path / "seed7.bam"
    output_b = tmp_path / "seed99.bam"
    reshape_bam(source_bam, reference_fasta, output_a, probabilities, seed=7)
    reshape_bam(source_bam, reference_fasta, output_b, probabilities, seed=99)

    def _qname_set(bam_path: Path) -> set[str]:
        with pysam.AlignmentFile(str(bam_path), "rb", reference_filename=str(reference_fasta)) as bam:
            return {read.query_name for read in bam.fetch(until_eof=True)}

    assert _qname_set(output_a) != _qname_set(output_b)


# --- histogram CSV round-trip -----------------------------------------------


def test_histogram_csv_round_trip(tmp_path: Path) -> None:
    histogram = Counter({100: 3, 150: 7, 500: 1})
    csv_path = write_histogram_csv(histogram, tmp_path / "hist.csv")

    loaded = read_histogram_csv(csv_path)

    assert loaded == histogram


def test_reshaping_config_accepts_precomputed_target_histogram(
    tmp_path: Path, reference_fasta: Path, source_bam: Path, target_bam: Path
) -> None:
    target_histogram = compute_insert_size_histogram(target_bam, reference_fasta)
    histogram_path = write_histogram_csv(target_histogram, tmp_path / "target_hist.csv")

    config = ReshapingConfig.from_sources(
        source_bam=source_bam,
        reference_fasta=reference_fasta,
        output_dir=tmp_path / "out2",
        target_histogram_path=histogram_path,
        ks_threshold=0.15,
        seed=42,
    )
    reshaper = FragmentReshaper(config)

    assert reshaper.load_target_histogram() == target_histogram
    result = reshaper.run()
    assert result["verification"].passed is True


def test_reshaping_config_rejects_both_or_neither_target_source() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        ReshapingConfig.from_sources(
            source_bam=Path("s.bam"), reference_fasta=Path("ref.fa"), output_dir=Path("out"),
        )
    with pytest.raises(ValueError, match="exactly one"):
        ReshapingConfig.from_sources(
            source_bam=Path("s.bam"),
            reference_fasta=Path("ref.fa"),
            output_dir=Path("out"),
            target_bam=Path("t.bam"),
            target_histogram_path=Path("t.csv"),
        )


# --- negative control: reshaped tumor vs. held-out cfDNA-like reads --------


def test_reshaped_tumor_is_not_trivially_separable_from_held_out_cfdna_reads(
    tmp_path: Path, reference_fasta: Path, source_bam: Path, target_bam: Path
) -> None:
    """The ticket's negative control, in its stated fallback form: "at minimum,
    that the KS/distance check passes on synthetic data proving the reshaping
    works". Here the check runs against a *held-out* cfDNA-like sample (a
    different seed from the same generating distribution, standing in for
    real_background's independent background reads) rather than against the
    exact histogram reshaping optimised against -- so passing is not
    tautological."""
    config = ReshapingConfig.from_sources(
        source_bam=source_bam,
        reference_fasta=reference_fasta,
        output_dir=tmp_path / "out",
        target_bam=target_bam,
        seed=42,
        sample_name="matched_normal_check",
    )
    reshaper = FragmentReshaper(config)
    result = reshaper.run(raise_on_failure=False)
    reshaped_histogram = read_histogram_csv(config.output_dir / f"{config.sample_name}.reshaped_histogram.csv")

    held_out_bam = _build_bam(tmp_path, "held_out_cfdna", _cfdna_like_tlens(4000, seed=999))
    held_out_histogram = compute_insert_size_histogram(held_out_bam, reference_fasta)

    held_out_check = verify_reshaping(reshaped_histogram, held_out_histogram, threshold=0.15)
    assert held_out_check.passed is True

    # Contrast: the *unreshaped* HCC1395-like source is trivially separable
    # from the same held-out cfDNA sample -- this is the failure mode the
    # ticket exists to fix, reproduced here as a sanity check on the test
    # itself.
    unreshaped_source_histogram = compute_insert_size_histogram(source_bam, reference_fasta)
    unreshaped_check = verify_reshaping(unreshaped_source_histogram, held_out_histogram, threshold=0.15)
    assert unreshaped_check.passed is False
    assert result["verification"].passed is True


# --- per-contig parallelism (histogram + retention passes) ----------------


MULTI_CONTIGS = ("chr1", "chr2", "chr3")


def _build_multi_contig_bam(
    tmp_path: Path,
    name: str,
    contig_tlens: dict[str, list[int]],
    cross_contig_pair: tuple[str, str] | None = None,
) -> Path:
    """A coordinate-sorted, multi-contig synthetic BAM: `contig_tlens` maps
    contig name -> list of normal, same-contig insert sizes (built with
    `_add_pair`, exactly like the single-contig fixtures above).

    `cross_contig_pair`, if given, is `(name, (contig_a, contig_b))` for one
    extra pair whose two mates are placed on *different* contigs -- read1 on
    `contig_a`, read2 on `contig_b` -- with `template_length=0` (as real
    aligners report for genuinely discordant pairs) and the proper-pair flag
    left *set* regardless, so `_passes_pair_filter`'s own flag check would
    otherwise let it through: this exercises `tlen == 0` as the actual
    guard that keeps a cross-contig pair out of both the histogram and the
    reshaped output, and does so identically whether that pair's two mates
    are fetched by the same per-contig worker or two different ones.
    """
    unsorted_path = tmp_path / f"{name}.unsorted.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": contig, "LN": CONTIG_LENGTH} for contig in contig_tlens],
    }
    contig_index = {contig: i for i, contig in enumerate(contig_tlens)}
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=header) as bam:
        for contig, tlens in contig_tlens.items():
            position = 1000
            for i, tlen in enumerate(tlens):
                _add_pair(bam, f"{name}_{contig}_{i}", position, tlen, reference_id=contig_index[contig])
                position += tlen + 20
        if cross_contig_pair is not None:
            pair_name, (contig_a, contig_b) = cross_contig_pair
            read1 = pysam.AlignedSegment(bam.header)
            read1.query_name = pair_name
            read1.query_sequence = "A" * READ_LENGTH
            read1.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
            read1.flag = PAIRED | PROPER_PAIR | READ1 | MATE_REVERSE
            read1.reference_id = contig_index[contig_a]
            read1.reference_start = 5000
            read1.mapping_quality = 60
            read1.cigar = [(0, READ_LENGTH)]
            read1.next_reference_id = contig_index[contig_b]
            read1.next_reference_start = 6000
            read1.template_length = 0
            read1.set_tag("MQ", 60)
            bam.write(read1)

            read2 = pysam.AlignedSegment(bam.header)
            read2.query_name = pair_name
            read2.query_sequence = "A" * READ_LENGTH
            read2.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
            read2.flag = PAIRED | PROPER_PAIR | READ2 | REVERSE
            read2.reference_id = contig_index[contig_b]
            read2.reference_start = 6000
            read2.mapping_quality = 60
            read2.cigar = [(0, READ_LENGTH)]
            read2.next_reference_id = contig_index[contig_a]
            read2.next_reference_start = 5000
            read2.template_length = 0
            read2.set_tag("MQ", 60)
            bam.write(read2)

    sorted_path = tmp_path / f"{name}.bam"
    pysam.sort("-o", str(sorted_path), str(unsorted_path))
    pysam.index(str(sorted_path))
    return sorted_path


def _multi_contig_reference_fasta(tmp_path: Path) -> Path:
    fasta_path = tmp_path / "multi_ref.fa"
    sequence = "ACGT" * (CONTIG_LENGTH // 4 + 1)
    with fasta_path.open("w") as handle:
        for contig in MULTI_CONTIGS:
            handle.write(f">{contig}\n")
            for i in range(0, CONTIG_LENGTH, 70):
                handle.write(sequence[i:i + 70] + "\n")
    pysam.faidx(str(fasta_path))
    return fasta_path


def _multi_contig_tlens() -> dict[str, list[int]]:
    return {
        "chr1": _hcc1395_like_tlens(600, seed=11),
        "chr2": _hcc1395_like_tlens(600, seed=12),
        "chr3": _hcc1395_like_tlens(600, seed=13),
    }


def test_parallel_histogram_matches_sequential_on_multi_contig_bam(tmp_path: Path) -> None:
    reference_fasta = _multi_contig_reference_fasta(tmp_path)
    bam_path = _build_multi_contig_bam(tmp_path, "multi_hist", _multi_contig_tlens())

    sequential = compute_insert_size_histogram(bam_path, reference_fasta, workers=1)
    parallel = compute_insert_size_histogram(bam_path, reference_fasta, workers=4)

    assert sum(sequential.values()) > 0
    assert parallel == sequential


def test_parallel_reshape_bam_matches_sequential_retained_reads(tmp_path: Path) -> None:
    reference_fasta = _multi_contig_reference_fasta(tmp_path)
    bam_path = _build_multi_contig_bam(tmp_path, "multi_reshape", _multi_contig_tlens())

    source_histogram = compute_insert_size_histogram(bam_path, reference_fasta, workers=1)
    # A synthetic, deliberately non-trivial target so acceptance probabilities
    # differ across the tlen range (not just "keep everything").
    target_histogram = Counter({tlen: max(1, 500 - tlen) for tlen in range(100, 501, 10)})
    probabilities = acceptance_probabilities(source_histogram, target_histogram)

    sequential_bam = tmp_path / "sequential.reshaped.bam"
    parallel_bam = tmp_path / "parallel.reshaped.bam"
    reshape_bam(bam_path, reference_fasta, sequential_bam, probabilities, seed=42, workers=1)
    reshape_bam(bam_path, reference_fasta, parallel_bam, probabilities, seed=42, workers=4)
    pysam.index(str(sequential_bam))
    pysam.index(str(parallel_bam))

    def _read_names(path: Path) -> list[str]:
        with pysam.AlignmentFile(str(path), "rb", reference_filename=str(reference_fasta)) as bam:
            return sorted(read.query_name for read in bam.fetch(until_eof=True))

    sequential_names = _read_names(sequential_bam)
    parallel_names = _read_names(parallel_bam)

    assert len(sequential_names) > 0
    assert parallel_names == sequential_names

    sequential_histogram = compute_insert_size_histogram(sequential_bam, reference_fasta, workers=1)
    parallel_histogram = compute_insert_size_histogram(parallel_bam, reference_fasta, workers=1)
    assert parallel_histogram == sequential_histogram


def test_cross_contig_mate_pair_excluded_identically_regardless_of_workers(tmp_path: Path) -> None:
    """A pair whose two mates land on different contigs (and are therefore
    fetched by two different per-contig workers under `workers>1`) must be
    excluded from both the histogram and the reshaped output the same way
    it would be sequentially -- never half-written (one mate kept, the other
    dropped) and never double-counted."""
    reference_fasta = _multi_contig_reference_fasta(tmp_path)
    bam_path = _build_multi_contig_bam(
        tmp_path,
        "cross_contig",
        _multi_contig_tlens(),
        cross_contig_pair=("discordant_pair", ("chr1", "chr3")),
    )

    sequential_histogram = compute_insert_size_histogram(bam_path, reference_fasta, workers=1)
    parallel_histogram = compute_insert_size_histogram(bam_path, reference_fasta, workers=4)
    assert parallel_histogram == sequential_histogram

    probabilities = {tlen: 1.0 for tlen in range(50, 501)}
    sequential_bam = tmp_path / "cross_sequential.bam"
    parallel_bam = tmp_path / "cross_parallel.bam"
    reshape_bam(bam_path, reference_fasta, sequential_bam, probabilities, seed=7, workers=1)
    reshape_bam(bam_path, reference_fasta, parallel_bam, probabilities, seed=7, workers=4)

    def _read_names(path: Path) -> list[str]:
        with pysam.AlignmentFile(str(path), "rb", reference_filename=str(reference_fasta)) as bam:
            return sorted(read.query_name for read in bam.fetch(until_eof=True))

    sequential_names = _read_names(sequential_bam)
    parallel_names = _read_names(parallel_bam)
    assert parallel_names == sequential_names
    # The discordant, tlen=0 cross-contig pair is excluded everywhere (both
    # its mates), not partially retained.
    assert "discordant_pair" not in sequential_names
    assert "discordant_pair" not in parallel_names
