"""Synthetic-BAM tests for the post-hoc reshaping pilot
(docs/todo/todo/reshaping-pilot-post-hoc-filter.md).

A synthetic "mixture" BAM carries two @RG groups, `TUMOR` (SM=HCC1395_WGS_NV_3,
a broad HCC1395-like insert-size spread) and `BACKGROUND` (SM=EGAFxxxx, a
tight cfDNA-like spread), mirroring what `admixture.py`'s merged/dedup BAMs
actually look like. Tests check: only tumor-RG pairs are thinned by the
insert-size acceptance table; both mates always agree; the result is
deterministic under a fixed seed; and the background-thinning step hits the
batch's common total-pair target within one pair per 10^4 (in fact exactly,
since `compute_background_keep_threshold` is an exact top-k cutoff, not an
independent per-pair probability).
"""

import importlib.util
import random
from collections import Counter
from pathlib import Path

import pysam
import pytest

from flash_fraq.fragment_reshaping import (
    acceptance_probabilities,
    compute_background_keep_threshold,
    count_mixture_reshape_pairs,
    mean_acceptance,
    reshape_mixture_bam,
    tumor_rg_ids,
    write_acceptance_table_csv,
    read_acceptance_table_csv,
)

SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "reshape_admixture_mixture.py"
_spec = importlib.util.spec_from_file_location("reshape_admixture_mixture", SCRIPT_PATH)
reshape_admixture_mixture = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(reshape_admixture_mixture)

CONTIG = "chr1"
CONTIG_LENGTH = 800_000
READ_LENGTH = 50

PAIRED = 0x1
PROPER_PAIR = 0x2
READ1 = 0x40
READ2 = 0x80
MATE_REVERSE = 0x20
REVERSE = 0x10

TUMOR_SM = "HCC1395_WGS_NV_3"
TUMOR_RG = "SRR_TUMOR_1"
BACKGROUND_SM = "EGAF00000001"
BACKGROUND_RG = "EGAF00000001"


def _add_pair(
    bam: pysam.AlignmentFile,
    name: str,
    start: int,
    tlen: int,
    rg: str,
    mapq1: int = 60,
    mapq2: int = 60,
) -> None:
    mate_start = start + tlen - READ_LENGTH

    read1 = pysam.AlignedSegment(bam.header)
    read1.query_name = name
    read1.query_sequence = "A" * READ_LENGTH
    read1.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
    read1.flag = PAIRED | PROPER_PAIR | READ1 | MATE_REVERSE
    read1.reference_id = 0
    read1.reference_start = start
    read1.mapping_quality = mapq1
    read1.cigar = [(0, READ_LENGTH)]
    read1.next_reference_id = 0
    read1.next_reference_start = mate_start
    read1.template_length = tlen
    read1.set_tag("RG", rg)
    read1.set_tag("MQ", mapq2)
    bam.write(read1)

    read2 = pysam.AlignedSegment(bam.header)
    read2.query_name = name
    read2.query_sequence = "A" * READ_LENGTH
    read2.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
    read2.flag = PAIRED | PROPER_PAIR | READ2 | REVERSE
    read2.reference_id = 0
    read2.reference_start = mate_start
    read2.mapping_quality = mapq2
    read2.cigar = [(0, READ_LENGTH)]
    read2.next_reference_id = 0
    read2.next_reference_start = start
    read2.template_length = -tlen
    read2.set_tag("RG", rg)
    read2.set_tag("MQ", mapq1)
    bam.write(read2)


def _header() -> dict:
    return {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": CONTIG, "LN": CONTIG_LENGTH}],
        "RG": [
            {"ID": TUMOR_RG, "SM": TUMOR_SM, "PL": "ILLUMINA"},
            {"ID": BACKGROUND_RG, "SM": BACKGROUND_SM, "PL": "HiSeq2500"},
        ],
    }


def _hcc1395_like_tlens(n: int, seed: int) -> list[int]:
    rng = random.Random(seed)
    return [rng.randint(100, 500) for _ in range(n)]


def _cfdna_like_tlens(n: int, seed: int) -> list[int]:
    rng = random.Random(seed)
    return [max(80, min(250, round(rng.gauss(166, 15)))) for _ in range(n)]


def _build_mixture_bam(
    tmp_path: Path, name: str, tumor_tlens: list[int], background_tlens: list[int]
) -> Path:
    unsorted_path = tmp_path / f"{name}.unsorted.bam"
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=_header()) as bam:
        position = 1000
        for i, tlen in enumerate(tumor_tlens):
            _add_pair(bam, f"tumor_{name}_{i}", position, tlen, TUMOR_RG)
            position += tlen + 20
        for i, tlen in enumerate(background_tlens):
            _add_pair(bam, f"bg_{name}_{i}", position, tlen, BACKGROUND_RG)
            position += tlen + 20
    sorted_path = tmp_path / f"{name}.dedup.bam"
    pysam.sort("-o", str(sorted_path), str(unsorted_path))
    pysam.index(str(sorted_path))
    return sorted_path


@pytest.fixture(scope="module")
def reference_fasta(tmp_path_factory) -> Path:
    tmp_path = tmp_path_factory.mktemp("ref")
    fasta_path = tmp_path / "ref.fa"
    sequence = "ACGT" * (CONTIG_LENGTH // 4 + 1)
    with fasta_path.open("w") as handle:
        handle.write(f">{CONTIG}\n")
        for i in range(0, CONTIG_LENGTH, 70):
            handle.write(sequence[i : i + 70] + "\n")
    pysam.faidx(str(fasta_path))
    return fasta_path


@pytest.fixture(scope="module")
def acceptance_table() -> dict[int, float]:
    source_counts = Counter(_hcc1395_like_tlens(4000, seed=11))
    target_counts = Counter(_cfdna_like_tlens(4000, seed=12))
    return acceptance_probabilities(source_counts, target_counts)


# --- tumor_rg_ids -------------------------------------------------------


def test_tumor_rg_ids_selects_only_the_tumor_sm(tmp_path: Path, reference_fasta: Path) -> None:
    bam_path = _build_mixture_bam(tmp_path, "mix", [200], [166])
    with pysam.AlignmentFile(str(bam_path), "rb", reference_filename=str(reference_fasta)) as bam:
        rg_ids = tumor_rg_ids(bam.header, TUMOR_SM)
    assert rg_ids == {TUMOR_RG}


# --- counting / reshaping: only tumor-RG pairs thinned -------------------


def test_only_tumor_rg_pairs_are_thinned(
    tmp_path: Path, reference_fasta: Path, acceptance_table: dict[int, float]
) -> None:
    tumor_tlens = _hcc1395_like_tlens(2000, seed=21)
    background_tlens = _cfdna_like_tlens(2000, seed=22)
    bam_path = _build_mixture_bam(tmp_path, "mix_a", tumor_tlens, background_tlens)

    counts = count_mixture_reshape_pairs(bam_path, reference_fasta, {TUMOR_RG}, acceptance_table, seed=7)

    assert counts.tumor_original == 2000
    assert counts.background_original == 2000
    # Some tumor pairs are actually thinned (acceptance table has both 0 and <1 entries).
    assert 0 < counts.tumor_kept < counts.tumor_original

    # Reshaping with background fully retained (threshold above every draw)
    # must keep every background pair and thin only the tumor RG.
    output_bam = tmp_path / "out" / "mix_a.reshaped.dedup.bam"
    _, final_counts = reshape_mixture_bam(
        bam_path, reference_fasta, output_bam, {TUMOR_RG}, acceptance_table, background_keep_threshold=1.0, seed=7
    )
    assert final_counts.background_kept == counts.background_original
    assert final_counts.tumor_kept == counts.tumor_kept

    with pysam.AlignmentFile(str(output_bam), "rb", reference_filename=str(reference_fasta)) as bam:
        kept_rgs = {read.get_tag("RG") for read in bam.fetch(until_eof=True)}
    assert kept_rgs == {TUMOR_RG, BACKGROUND_RG}


# --- MAPQ filter: tumor only, background untouched ------------------------


def test_low_mapq_tumor_pairs_are_excluded_but_background_is_not(
    tmp_path: Path, reference_fasta: Path
) -> None:
    unsorted_path = tmp_path / "mix_mapq.unsorted.bam"
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=_header()) as bam:
        _add_pair(bam, "tumor_high", 1000, 200, TUMOR_RG, mapq1=60, mapq2=60)
        _add_pair(bam, "tumor_low", 2000, 200, TUMOR_RG, mapq1=10, mapq2=10)
        _add_pair(bam, "bg_high", 3000, 166, BACKGROUND_RG, mapq1=60, mapq2=60)
        # A low-MAPQ background pair must survive untouched -- background
        # reads are never MAPQ-filtered here, only tumor-RG reads are.
        _add_pair(bam, "bg_low", 4000, 166, BACKGROUND_RG, mapq1=5, mapq2=5)
    bam_path = tmp_path / "mix_mapq.dedup.bam"
    pysam.sort("-o", str(bam_path), str(unsorted_path))
    pysam.index(str(bam_path))

    acceptance_table = {200: 1.0, 166: 1.0}

    counts = count_mixture_reshape_pairs(
        bam_path, reference_fasta, {TUMOR_RG}, acceptance_table, seed=1, min_mapq=30
    )
    assert counts.tumor_original == 2
    assert counts.tumor_kept == 1  # only tumor_high
    assert counts.background_original == 2  # both background pairs counted, untouched here

    output_bam = tmp_path / "out" / "mix_mapq.reshaped.dedup.bam"
    _, final_counts = reshape_mixture_bam(
        bam_path, reference_fasta, output_bam, {TUMOR_RG}, acceptance_table,
        background_keep_threshold=1.0, seed=1, min_mapq=30,
    )
    assert final_counts.tumor_kept == 1
    assert final_counts.background_kept == 2  # bg_low kept despite MAPQ 5

    with pysam.AlignmentFile(str(output_bam), "rb", reference_filename=str(reference_fasta)) as bam:
        kept_names = {read.query_name for read in bam.fetch(until_eof=True)}
    assert kept_names == {"tumor_high", "bg_high", "bg_low"}


def test_min_mapq_threshold_is_configurable_for_mixtures(tmp_path: Path, reference_fasta: Path) -> None:
    unsorted_path = tmp_path / "mix_mapq_cfg.unsorted.bam"
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=_header()) as bam:
        _add_pair(bam, "tumor_40", 1000, 200, TUMOR_RG, mapq1=40, mapq2=40)
    bam_path = tmp_path / "mix_mapq_cfg.dedup.bam"
    pysam.sort("-o", str(bam_path), str(unsorted_path))
    pysam.index(bam_path.as_posix())

    acceptance_table = {200: 1.0}

    kept_at_30 = count_mixture_reshape_pairs(bam_path, reference_fasta, {TUMOR_RG}, acceptance_table, seed=1, min_mapq=30)
    kept_at_50 = count_mixture_reshape_pairs(bam_path, reference_fasta, {TUMOR_RG}, acceptance_table, seed=1, min_mapq=50)
    assert kept_at_30.tumor_kept == 1
    assert kept_at_50.tumor_kept == 0


# --- mates always agree --------------------------------------------------


def test_mates_always_agree(tmp_path: Path, reference_fasta: Path, acceptance_table: dict[int, float]) -> None:
    tumor_tlens = _hcc1395_like_tlens(1500, seed=31)
    background_tlens = _cfdna_like_tlens(1500, seed=32)
    bam_path = _build_mixture_bam(tmp_path, "mix_b", tumor_tlens, background_tlens)
    threshold = compute_background_keep_threshold(bam_path, reference_fasta, {TUMOR_RG}, background_keep_count=750, seed=9)

    output_bam = tmp_path / "out" / "mix_b.reshaped.dedup.bam"
    reshape_mixture_bam(bam_path, reference_fasta, output_bam, {TUMOR_RG}, acceptance_table, threshold, seed=9)

    with pysam.AlignmentFile(str(output_bam), "rb", reference_filename=str(reference_fasta)) as bam:
        names = [read.query_name for read in bam.fetch(until_eof=True)]
    assert len(names) > 0
    assert all(count == 2 for count in Counter(names).values())


# --- determinism -----------------------------------------------------------


def test_reshape_mixture_bam_is_deterministic_under_fixed_seed(
    tmp_path: Path, reference_fasta: Path, acceptance_table: dict[int, float]
) -> None:
    tumor_tlens = _hcc1395_like_tlens(1000, seed=41)
    background_tlens = _cfdna_like_tlens(1000, seed=42)
    bam_path = _build_mixture_bam(tmp_path, "mix_c", tumor_tlens, background_tlens)
    threshold = compute_background_keep_threshold(bam_path, reference_fasta, {TUMOR_RG}, background_keep_count=400, seed=5)

    out_a = tmp_path / "a.bam"
    out_b = tmp_path / "b.bam"
    reshape_mixture_bam(bam_path, reference_fasta, out_a, {TUMOR_RG}, acceptance_table, threshold, seed=5)
    reshape_mixture_bam(bam_path, reference_fasta, out_b, {TUMOR_RG}, acceptance_table, threshold, seed=5)

    def _qnames(path: Path) -> list[str]:
        with pysam.AlignmentFile(str(path), "rb", reference_filename=str(reference_fasta)) as bam:
            return sorted(read.query_name for read in bam.fetch(until_eof=True))

    assert _qnames(out_a) == _qnames(out_b)


def test_reshape_mixture_bam_seed_changes_retained_background_set(
    tmp_path: Path, reference_fasta: Path, acceptance_table: dict[int, float]
) -> None:
    tumor_tlens = _hcc1395_like_tlens(1000, seed=51)
    background_tlens = _cfdna_like_tlens(1000, seed=52)
    bam_path = _build_mixture_bam(tmp_path, "mix_d", tumor_tlens, background_tlens)

    threshold_a = compute_background_keep_threshold(bam_path, reference_fasta, {TUMOR_RG}, 400, seed=1)
    threshold_b = compute_background_keep_threshold(bam_path, reference_fasta, {TUMOR_RG}, 400, seed=2)
    out_a = tmp_path / "seed1.bam"
    out_b = tmp_path / "seed2.bam"
    reshape_mixture_bam(bam_path, reference_fasta, out_a, {TUMOR_RG}, acceptance_table, threshold_a, seed=1)
    reshape_mixture_bam(bam_path, reference_fasta, out_b, {TUMOR_RG}, acceptance_table, threshold_b, seed=2)

    def _bg_qnames(path: Path) -> set[str]:
        with pysam.AlignmentFile(str(path), "rb", reference_filename=str(reference_fasta)) as bam:
            return {read.query_name for read in bam.fetch(until_eof=True) if read.get_tag("RG") == BACKGROUND_RG}

    assert _bg_qnames(out_a) != _bg_qnames(out_b)


# --- exact background quota -------------------------------------------------


def test_background_keep_threshold_hits_exact_target(
    tmp_path: Path, reference_fasta: Path, acceptance_table: dict[int, float]
) -> None:
    tumor_tlens = _hcc1395_like_tlens(3000, seed=61)
    background_tlens = _cfdna_like_tlens(7000, seed=62)
    bam_path = _build_mixture_bam(tmp_path, "mix_e", tumor_tlens, background_tlens)

    for keep_count in (0, 1, 500, 3333, 7000, 8000):
        threshold = compute_background_keep_threshold(bam_path, reference_fasta, {TUMOR_RG}, keep_count, seed=13)
        output_bam = tmp_path / f"out_{keep_count}.bam"
        _, final_counts = reshape_mixture_bam(
            bam_path, reference_fasta, output_bam, {TUMOR_RG}, acceptance_table, threshold, seed=13
        )
        expected = min(keep_count, 7000)
        assert final_counts.background_kept == expected


# --- one-pair-per-10^4 batch normalisation, via the CLI script ------------


def test_script_normalises_total_pairs_to_a_common_target(tmp_path: Path, reference_fasta: Path) -> None:
    """Three mixtures at different (synthetic) tumor fractions -- more tumor
    reads means more thinning by the acceptance table, so the raw totals
    differ. After the script's background-thinning pass, every mixture's
    total kept-pair count must match the common target within 1 pair per
    10^4 (here: exactly, since the cutoff is an exact top-k selection)."""
    source_counts = Counter(_hcc1395_like_tlens(20000, seed=100))
    target_counts = Counter(_cfdna_like_tlens(20000, seed=101))
    table = acceptance_probabilities(source_counts, target_counts)
    table_path = write_acceptance_table_csv(table, tmp_path / "acceptance.csv")

    mixtures_dir = tmp_path / "mixtures"
    mixtures_dir.mkdir()
    n_tumor_per_mixture = {"low": 500, "mid": 2000, "high": 5000}
    n_background = 5000
    bam_paths = []
    for i, (label, n_tumor) in enumerate(n_tumor_per_mixture.items()):
        tumor_tlens = _hcc1395_like_tlens(n_tumor, seed=200 + i)
        background_tlens = _cfdna_like_tlens(n_background, seed=300 + i)
        bam_paths.append(_build_mixture_bam(mixtures_dir, f"mix_{label}", tumor_tlens, background_tlens))

    output_dir = tmp_path / "reshaped"
    argv = [
        "reshape_admixture_mixture.py",
        "--tumor-sample-name", TUMOR_SM,
        "--acceptance-table", str(table_path),
        "--reference-fasta", str(reference_fasta),
        "--output-dir", str(output_dir),
        "--seed", "77",
    ]
    for bam_path in bam_paths:
        argv.extend(["--mixture-bam", str(bam_path)])

    import sys

    old_argv = sys.argv
    sys.argv = argv
    try:
        reshape_admixture_mixture.main()
    finally:
        sys.argv = old_argv

    import csv as csv_module

    with (output_dir / "reshape_admixture_mixture_summary.csv").open() as handle:
        rows = list(csv_module.DictReader(handle))

    assert len(rows) == 3
    common_target = int(rows[0]["common_target_total"])
    for row in rows:
        total_kept = int(row["total_kept"])
        # "within one pair per 10^4" -- common_target_total is itself in the
        # thousands here, so a tolerance of max(1, target * 1e-4) covers both
        # the ticket's relative bound and small-N exactness.
        tolerance = max(1, round(common_target * 1e-4))
        assert abs(total_kept - common_target) <= tolerance
        assert Path(row["output_bam"]).exists()
        # Never overwrites the original dedup.bam.
        assert row["output_bam"] != row["input_bam"]
        assert Path(row["output_bam"]).name.endswith(".reshaped.dedup.bam")


def test_cli_min_mapq_flag_drops_low_mapq_tumor_reads(tmp_path: Path, reference_fasta: Path) -> None:
    source_counts = Counter(_hcc1395_like_tlens(500, seed=500))
    target_counts = Counter(_cfdna_like_tlens(500, seed=501))
    table = acceptance_probabilities(source_counts, target_counts)
    table_path = write_acceptance_table_csv(table, tmp_path / "acceptance.csv")

    mixtures_dir = tmp_path / "mixtures"
    mixtures_dir.mkdir()
    unsorted_path = mixtures_dir / "mix_cli.unsorted.bam"
    with pysam.AlignmentFile(str(unsorted_path), "wb", header=_header()) as bam:
        _add_pair(bam, "tumor_high", 1000, 166, TUMOR_RG, mapq1=60, mapq2=60)
        _add_pair(bam, "tumor_low", 2000, 166, TUMOR_RG, mapq1=10, mapq2=10)
    bam_path = mixtures_dir / "mix_cli.dedup.bam"
    pysam.sort("-o", str(bam_path), str(unsorted_path))
    pysam.index(str(bam_path))

    output_dir = tmp_path / "reshaped_cli"
    argv = [
        "reshape_admixture_mixture.py",
        "--tumor-sample-name", TUMOR_SM,
        "--acceptance-table", str(table_path),
        "--reference-fasta", str(reference_fasta),
        "--output-dir", str(output_dir),
        "--seed", "77",
        "--min-mapq", "30",
        "--mixture-bam", str(bam_path),
    ]

    import sys

    old_argv = sys.argv
    sys.argv = argv
    try:
        reshape_admixture_mixture.main()
    finally:
        sys.argv = old_argv

    output_bam = output_dir / "mix_cli.reshaped.dedup.bam"
    with pysam.AlignmentFile(str(output_bam), "rb", reference_filename=str(reference_fasta)) as bam:
        kept_names = {read.query_name for read in bam.fetch(until_eof=True)}
    assert "tumor_low" not in kept_names


def test_output_dir_must_not_collide_with_source_dir(tmp_path: Path, reference_fasta: Path) -> None:
    source_counts = Counter(_hcc1395_like_tlens(200, seed=400))
    target_counts = Counter(_cfdna_like_tlens(200, seed=401))
    table = acceptance_probabilities(source_counts, target_counts)
    table_path = write_acceptance_table_csv(table, tmp_path / "acceptance.csv")

    mixtures_dir = tmp_path / "mixtures"
    mixtures_dir.mkdir()
    bam_path = _build_mixture_bam(
        mixtures_dir, "mix_only", _hcc1395_like_tlens(200, seed=402), _cfdna_like_tlens(200, seed=403)
    )

    import sys

    old_argv = sys.argv
    sys.argv = [
        "reshape_admixture_mixture.py",
        "--tumor-sample-name", TUMOR_SM,
        "--acceptance-table", str(table_path),
        "--reference-fasta", str(reference_fasta),
        "--output-dir", str(mixtures_dir),
        "--mixture-bam", str(bam_path),
    ]
    try:
        with pytest.raises(ValueError, match="own directory"):
            reshape_admixture_mixture.main()
    finally:
        sys.argv = old_argv


# --- acceptance-table CSV round trip / mean_acceptance ----------------------


def test_acceptance_table_csv_round_trip(tmp_path: Path, acceptance_table: dict[int, float]) -> None:
    csv_path = write_acceptance_table_csv(acceptance_table, tmp_path / "table.csv")
    loaded = read_acceptance_table_csv(csv_path)
    assert loaded.keys() == acceptance_table.keys()
    for insert_size, probability in acceptance_table.items():
        assert loaded[insert_size] == pytest.approx(probability)


def test_mean_acceptance_is_between_zero_and_one(acceptance_table: dict[int, float]) -> None:
    source_counts = Counter(_hcc1395_like_tlens(4000, seed=11))
    a_bar = mean_acceptance(source_counts, acceptance_table)
    assert 0.0 < a_bar <= 1.0
