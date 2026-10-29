"""In silico fragment-size reshaping (docs/todo/fragment-size-reshaping.md).

HCC1395/HCC1395BL is mechanically sheared cell-line WGS with a roughly
uniform/Gaussian insert-size distribution; real ctDNA carries the
short-fragment enrichment of *in vivo* nuclease digestion, which is this
project's entire feature set. Mixing HCC1395 into a cfDNA background
unreshaped would teach a "sheared vs. digested library" classifier instead of
real tumor-fraction signal (see the ticket's "Why this exists"). The ticket
names two mitigations and the project owner picked the first: reshape
HCC1395's fragment sizes to look like real ctDNA *before* mixing.

This is "in silico fragmentation" by insert-size-dependent read-pair
retention (rejection sampling), not re-fragmentation: it can only ever
*downsample* a source's own existing fragments toward a target shape, never
invent fragments a source lacks entirely. It operates on the already-aligned,
coordinate-sorted BAM/CRAM (the ticket's explicitly preferred design fork),
so it needs no HCC1395 realignment and no re-derived FASTQ.

Three stages:

1. Insert-size histograms, both target (drawn from a config-supplied
   background BAM/CRAM's own observed fragment sizes -- never a hardcoded
   shape) and source (the tumor BAM being reshaped). `compute_insert_size_histogram`
   reuses `fragment_extraction`'s own properly-paired/non-duplicate read-pair
   predicate (`REQUIRE_FLAGS`/`EXCLUDE_FLAGS`/`READ1_FLAG`) rather than
   re-deriving that logic.
   Both the histogram and every later thinning step also require a minimum
   MAPQ (`DEFAULT_MIN_MAPQ`, both mates via `_meets_mapq_threshold`):
   `coverage_bias_report.json` found MAPQ<30 pairs accepted 6-18x more often
   than MAPQ>=30 ones under tlen-only rejection sampling, a repeat/
   mappability artefact. This is applied to the tumor source everywhere and
   to the background/target only when computing histograms (never when
   thinning an already-merged mixture's real background reads -- those are
   not this artefact's source).
2. `acceptance_probabilities`: one retention probability per insert-size bin,
   `target_density / (M * source_density)` where `M` is chosen so every
   probability is <= 1 -- standard rejection-sampling reweighting. A read
   pair's accept/reject draw is deterministic given `(seed, qname, insert
   size)`, so both mates of a pair always agree and a re-run with the same
   seed reproduces the same retained set exactly.
3. `verify_reshaping`: the reshaped population's insert-size histogram is
   compared against the target via a weighted Kolmogorov-Smirnov statistic
   (the max gap between the two empirical CDFs, computed directly from
   histograms rather than by materialising per-read samples -- the ticket's
   "KS test or histogram-distance" option, picked because it scales to a real
   BAM's fragment counts without ever holding a per-fragment list in memory).

`FragmentReshaper` drives all three end to end. It produces a reshaped
BAM that a later wiring step would hand to `admixture.py` as an alternative
tumor source -- `admixture.py`'s mixing logic itself is untouched here.
"""

from __future__ import annotations

import csv
import hashlib
import json
import tempfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import pysam

from flash_fraq.fragment_extraction import EXCLUDE_FLAGS, READ1_FLAG, REQUIRE_FLAGS

# ichorCNA-style ULP-WGS libraries and cfDNA fragments both live well inside
# this range; anything outside it is either a mapping artefact or a rare tail
# not worth histogram resolution over. Configurable per `ReshapingConfig`.
DEFAULT_MIN_INSERT_SIZE: int = 30
DEFAULT_MAX_INSERT_SIZE: int = 700

# `ks_statistic_from_histograms` output is a CDF-distance in [0, 1]; 0.15 is a
# deliberately loose default so a first real run fails loud on genuine mismatch
# rather than on synthetic-scale sampling noise. Override per invocation.
DEFAULT_KS_THRESHOLD: float = 0.15

# reshaping_pilot/coverage_check/coverage_bias_report.json found MAPQ<30 read
# pairs accepted 6-18x more often than MAPQ>=30 pairs -- a repeat/mappability
# artefact of rejection sampling on tlen alone, since low-MAPQ multi-mapping
# regions skew toward the insert sizes this module's acceptance table favours.
# 30 is the conventional "confidently uniquely mapped" cutoff (same one
# `ichorCNA`/most WGS QC pipelines use) and matches what the coverage-bias
# report itself was computed against.
DEFAULT_MIN_MAPQ: int = 30

HISTOGRAM_CSV_COLUMNS: tuple[str, ...] = ("insert_size", "count")
ACCEPTANCE_TABLE_CSV_COLUMNS: tuple[str, ...] = ("insert_size", "probability")

# The uniform, tlen-independent background-thinning draw
# (`reshape_admixture_mixture.py`'s `k`) reuses `_acceptance_draw`'s hash with
# this sentinel in place of a real insert size -- real `abs(template_length)`
# values are never negative, so this can never collide with a genuine tumor
# insert-size draw for the same qname/seed.
BACKGROUND_THIN_SENTINEL: int = -1


def _passes_pair_filter(read: pysam.AlignedSegment) -> bool:
    """Properly paired, not duplicate/secondary/supplementary -- both mates.

    Deliberately omits `fragment_extraction`'s `READ1_FLAG` check: that one is
    for counting each fragment once in a *histogram*, not for deciding which
    BAM records participate in reshaping (both mates of a retained pair must
    be written).
    """
    if read.flag & EXCLUDE_FLAGS:
        return False
    if not (read.flag & REQUIRE_FLAGS):
        return False
    return True


def _meets_mapq_threshold(read: pysam.AlignedSegment, min_mapq: int) -> bool:
    """Both-mates MAPQ gate for one read record.

    A single pass over a coordinate-sorted BAM sees each mate as its own
    record with no cross-mate state, yet both mates of a retained pair must
    reach the *same* accept/reject decision (`reshape_mixture_bam` writes
    each mate independently). This record's own `mapping_quality` is exactly
    that mate's MAPQ; the other mate's MAPQ is only available via the `MQ`
    tag, which `bwa mem` (this project's aligner) writes on every properly
    paired read. Using `min(own MAPQ, MQ tag)` is symmetric: read1's record
    computes `min(mapq_r1, mapq_r2)` and read2's record computes
    `min(mapq_r2, mapq_r1)` -- the same value -- so both mates independently
    arrive at the identical decision without needing to see each other.

    If `MQ` is absent (older/other aligners), this falls back to the
    record's own MAPQ alone; that is a best-effort approximation only (the
    two mates of such a pair could then disagree), documented here rather
    than treated as an error since `bwa mem` output is this project's only
    real input.
    """
    own_mapq = read.mapping_quality
    mate_mapq = read.get_tag("MQ") if read.has_tag("MQ") else own_mapq
    return min(own_mapq, mate_mapq) >= min_mapq


def _reference_names_with_reads(bam_path: Path, reference_fasta: Path) -> list[str]:
    """Contig names to hand out as per-worker fetch units for the parallel
    histogram/retention passes.

    Prefers the BAM index's own per-contig mapped/unmapped counts (cheap,
    already on disk since every input here is indexed for `pysam.fetch`) to
    skip contigs with zero reads entirely -- on hg19 that is most of the
    unplaced/random/decoy contigs. Falls back to the full header contig list
    if index stats are unavailable (e.g. a CRAM without a compatible index),
    since an empty `bam.fetch(contig)` on a real contig is cheap, just not
    free.
    """
    with pysam.AlignmentFile(str(bam_path), "rb", reference_filename=str(reference_fasta)) as bam:
        try:
            stats = bam.get_index_statistics()
            names = [s.contig for s in stats if (s.mapped + s.unmapped) > 0]
            if names:
                return names
        except (ValueError, AttributeError):
            pass
        return list(bam.references)


def _histogram_worker(
    bam_path: str,
    reference_fasta: str,
    contig: str,
    min_insert_size: int,
    max_insert_size: int,
    min_mapq: int,
) -> Counter:
    """One contig's worth of `compute_insert_size_histogram`'s sequential
    body, run in its own process. `_passes_pair_filter`'s `REQUIRE_FLAGS`
    (proper-pair, 0x2) already implies both mates are mapped on the same
    contig -- an aligner only sets it for such pairs -- so fetching by
    contig here can never split a qualifying pair across two workers or
    silently drop an unmapped/cross-contig one: those are excluded by the
    same filter regardless of whether the pass is sequential or per-contig.
    """
    histogram: Counter = Counter()
    with pysam.AlignmentFile(bam_path, "rb", reference_filename=reference_fasta) as bam:
        for read in bam.fetch(contig):
            if not _passes_pair_filter(read):
                continue
            if not (read.flag & READ1_FLAG):
                continue
            if not _meets_mapq_threshold(read, min_mapq):
                continue
            tlen = abs(read.template_length)
            if tlen == 0:
                continue
            if min_insert_size <= tlen <= max_insert_size:
                histogram[tlen] += 1
    return histogram


def compute_insert_size_histogram(
    bam_path: Path,
    reference_fasta: Path,
    min_insert_size: int = DEFAULT_MIN_INSERT_SIZE,
    max_insert_size: int = DEFAULT_MAX_INSERT_SIZE,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    workers: int = 1,
) -> Counter:
    """One count per fragment (read1 of each properly-paired, non-duplicate
    pair), keyed by `abs(template_length)`, restricted to `[min, max]`.

    `min_mapq` (see `_meets_mapq_threshold`) is applied uniformly to
    whichever BAM this is called against -- source (tumor) or target
    (background) -- so the acceptance table built from their ratio compares
    like-for-like populations; see `coverage_bias_report.json`'s finding that
    low-MAPQ reads are accepted at a very different rate than high-MAPQ ones
    under tlen-only rejection sampling.

    `reference_filename` is required for CRAM decoding on these files -- same
    reasoning as `fragment_extraction.FragmentExtractor.extract_midpoints`.

    `workers <= 1` (the default) runs the original single-pass,
    `fetch(until_eof=True)` loop unchanged -- byte-for-byte the same code
    path as before this parameter existed, so existing callers and outputs
    are untouched. `workers > 1` instead fans the same per-read predicate out
    per contig across a process pool and sums the per-contig `Counter`s; the
    per-pair accept/reject logic itself is untouched (this is a histogram
    tally, not a thinning decision) so the summed result is the same
    multiset of counts regardless of how many workers computed it, just
    order-independent-summed rather than accumulated in file order.
    """
    if workers <= 1:
        histogram: Counter = Counter()
        with pysam.AlignmentFile(str(bam_path), "rb", reference_filename=str(reference_fasta)) as bam:
            for read in bam.fetch(until_eof=True):
                if not _passes_pair_filter(read):
                    continue
                if not (read.flag & READ1_FLAG):
                    continue
                if not _meets_mapq_threshold(read, min_mapq):
                    continue
                tlen = abs(read.template_length)
                if tlen == 0:
                    continue
                if min_insert_size <= tlen <= max_insert_size:
                    histogram[tlen] += 1
        return histogram

    contigs = _reference_names_with_reads(bam_path, reference_fasta)
    histogram = Counter()
    if not contigs:
        return histogram
    with ProcessPoolExecutor(max_workers=min(workers, len(contigs))) as pool:
        futures = [
            pool.submit(
                _histogram_worker,
                str(bam_path),
                str(reference_fasta),
                contig,
                min_insert_size,
                max_insert_size,
                min_mapq,
            )
            for contig in contigs
        ]
        for future in futures:
            histogram.update(future.result())
    return histogram


def write_histogram_csv(histogram: Counter, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(HISTOGRAM_CSV_COLUMNS)
        for insert_size in sorted(histogram):
            writer.writerow([insert_size, histogram[insert_size]])
    return output_path


def read_histogram_csv(input_path: Path) -> Counter:
    histogram: Counter = Counter()
    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            histogram[int(row["insert_size"])] = int(row["count"])
    return histogram


def write_acceptance_table_csv(probabilities: dict[int, float], output_path: Path) -> Path:
    """Round-trip format for a precomputed `acceptance_probabilities` table,
    mirroring `write_histogram_csv` -- so an expensive real-data acceptance
    table (pooled-background target x HCC1395 source) can be computed once
    and reused across every mixture `reshape_admixture_mixture.py` touches."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(ACCEPTANCE_TABLE_CSV_COLUMNS)
        for insert_size in sorted(probabilities):
            writer.writerow([insert_size, probabilities[insert_size]])
    return output_path


def read_acceptance_table_csv(input_path: Path) -> dict[int, float]:
    probabilities: dict[int, float] = {}
    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            probabilities[int(row["insert_size"])] = float(row["probability"])
    return probabilities


def mean_acceptance(source_histogram: Counter, probabilities: dict[int, float]) -> float:
    """The source-weighted mean acceptance probability (script/ticket's
    "a-bar") -- the expected fraction of tumor pairs retained overall, as
    opposed to `probabilities`' one rate per insert size."""
    source_total = sum(source_histogram.values())
    if source_total == 0:
        raise ValueError("source histogram is empty; nothing to average acceptance over")
    weighted = sum(count * probabilities.get(insert_size, 0.0) for insert_size, count in source_histogram.items())
    return weighted / source_total


def tumor_rg_ids(bam_header: "pysam.AlignmentHeader | dict", tumor_sample_name: str) -> set[str]:
    """RG IDs in a merged mixture BAM's header whose `SM` tag equals
    `tumor_sample_name` -- `admixture.py`'s merged/dedup BAMs keep each
    source's own distinct `@RG` `SM` tag (see its `compute_empirical_tumor_fraction`
    docstring), so this is how `reshape_admixture_mixture.py` tells tumor reads
    from background reads apart in an already-merged BAM without needing a
    second, unmerged copy of either source."""
    header_dict = bam_header.to_dict() if hasattr(bam_header, "to_dict") else bam_header
    return {rg["ID"] for rg in header_dict.get("RG", []) if rg.get("SM") == tumor_sample_name}


def background_thin_draw(seed: int, qname: str) -> float:
    """A deterministic, tlen-independent draw in [0, 1) for background-pair
    thinning -- same construction as `_acceptance_draw`, keyed by the
    `BACKGROUND_THIN_SENTINEL` sentinel instead of a real insert size."""
    return _acceptance_draw(seed, qname, BACKGROUND_THIN_SENTINEL)


def acceptance_probabilities(source_histogram: Counter, target_histogram: Counter) -> dict[int, float]:
    """One retention probability per insert size present in the source.

    Standard rejection-sampling reweighting: `ratio(size) = target_density(size)
    / source_density(size)`, then every probability is scaled down by the
    largest ratio so none exceeds 1 (`samtools view -s` logic all over again --
    a distribution can only ever be thinned, never inflated). An insert size
    with zero target density (a length no real cfDNA fragment in the target
    background ever has) gets probability 0: reshaping can suppress a source's
    own sizes but can never manufacture a size it lacks.
    """
    source_total = sum(source_histogram.values())
    target_total = sum(target_histogram.values())
    if source_total == 0:
        raise ValueError("source histogram is empty; nothing to reshape")
    if target_total == 0:
        raise ValueError("target histogram is empty; nothing to reshape toward")

    ratios: dict[int, float] = {}
    for insert_size, source_count in source_histogram.items():
        source_density = source_count / source_total
        target_density = target_histogram.get(insert_size, 0) / target_total
        ratios[insert_size] = target_density / source_density if source_density > 0 else 0.0

    max_ratio = max(ratios.values(), default=0.0)
    if max_ratio <= 0:
        raise ValueError(
            "source and target insert-size histograms share no overlapping bins; "
            "reshaping cannot retain any read pairs"
        )
    return {insert_size: ratio / max_ratio for insert_size, ratio in ratios.items()}


def _acceptance_draw(seed: int, qname: str, insert_size: int) -> float:
    """A deterministic draw in [0, 1) from `(seed, qname, insert_size)`.

    Keyed by the read pair's own identity rather than by iteration order, so
    both mates of a pair -- which may appear anywhere relative to each other
    in a coordinate-sorted BAM -- independently compute the identical
    accept/reject decision with no cross-mate state to track, and a re-run
    with the same seed reproduces the same retained set exactly.
    """
    digest = hashlib.sha256(f"{seed}:{qname}:{insert_size}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") / float(1 << 64)


def _reshape_contig_worker(
    source_bam: str,
    reference_fasta: str,
    contig: str,
    acceptance_probabilities: dict[int, float],
    seed: int,
    min_insert_size: int,
    max_insert_size: int,
    min_mapq: int,
    output_path: str,
) -> str:
    """One contig's worth of `reshape_bam`'s sequential retention body,
    writing a standalone per-contig BAM (same full header as the source, so
    every per-contig file is independently valid and header-compatible for
    `pysam.cat`). Same reasoning as `_histogram_worker` for why per-contig
    fetch cannot split or drop a qualifying pair across workers. The
    accept/reject draw (`_acceptance_draw(seed, qname, tlen)`) depends only
    on `(seed, qname, tlen)`, never on worker count or iteration order, so
    this worker makes exactly the same decision for a given read that the
    sequential path would.
    """
    with pysam.AlignmentFile(source_bam, "rb", reference_filename=reference_fasta) as bam_in:
        with pysam.AlignmentFile(output_path, "wb", header=bam_in.header) as bam_out:
            for read in bam_in.fetch(contig):
                if not _passes_pair_filter(read):
                    continue
                if not _meets_mapq_threshold(read, min_mapq):
                    continue
                tlen = abs(read.template_length)
                if tlen == 0 or not (min_insert_size <= tlen <= max_insert_size):
                    continue
                probability = acceptance_probabilities.get(tlen, 0.0)
                if probability <= 0:
                    continue
                if probability < 1.0 and _acceptance_draw(seed, read.query_name, tlen) >= probability:
                    continue
                bam_out.write(read)
    return output_path


def reshape_bam(
    source_bam: Path,
    reference_fasta: Path,
    output_bam: Path,
    acceptance_probabilities: dict[int, float],
    seed: int,
    min_insert_size: int = DEFAULT_MIN_INSERT_SIZE,
    max_insert_size: int = DEFAULT_MAX_INSERT_SIZE,
    min_mapq: int = DEFAULT_MIN_MAPQ,
    workers: int = 1,
) -> Path:
    """Writes the retained read pairs to `output_bam`, coordinate-sorted.

    `source_bam` here is always the tumor source being reshaped (this
    function has no background/tumor distinction to make), so `min_mapq`
    (see `_meets_mapq_threshold`) is applied to every pair unconditionally --
    below-threshold tumor pairs are dropped before the tlen acceptance draw
    even runs.

    `workers <= 1` (the default): unchanged single pass, `fetch(until_eof=
    True)` yields records in file order (the source is already
    coordinate-sorted, per this module's post-alignment design), and a
    subset of a sorted sequence is still sorted, so no re-sort is needed
    before indexing. Byte-for-byte the same code path as before this
    parameter existed.

    `workers > 1`: the same per-read retention predicate runs per contig in
    a process pool, each worker writing its own per-contig BAM
    (`_reshape_contig_worker`), then the per-contig outputs are concatenated
    in header contig order via `pysam.cat` (not `samtools merge` -- inputs
    are disjoint contig ranges, not overlapping records needing a merge
    sort) and the combined file is indexed. Concatenating disjoint,
    individually-sorted contig ranges in header order reproduces the same
    global coordinate order as the sequential path, and the retained set of
    read pairs (and their per-read record bytes) is identical regardless of
    worker count since `_acceptance_draw` keys only on `(seed, qname,
    tlen)`.
    """
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    if workers <= 1:
        with pysam.AlignmentFile(str(source_bam), "rb", reference_filename=str(reference_fasta)) as bam_in:
            with pysam.AlignmentFile(str(output_bam), "wb", header=bam_in.header) as bam_out:
                for read in bam_in.fetch(until_eof=True):
                    if not _passes_pair_filter(read):
                        continue
                    if not _meets_mapq_threshold(read, min_mapq):
                        continue
                    tlen = abs(read.template_length)
                    if tlen == 0 or not (min_insert_size <= tlen <= max_insert_size):
                        continue
                    probability = acceptance_probabilities.get(tlen, 0.0)
                    if probability <= 0:
                        continue
                    if probability < 1.0 and _acceptance_draw(seed, read.query_name, tlen) >= probability:
                        continue
                    bam_out.write(read)
        return output_bam

    contigs = _reference_names_with_reads(source_bam, reference_fasta)
    if not contigs:
        with pysam.AlignmentFile(str(source_bam), "rb", reference_filename=str(reference_fasta)) as bam_in:
            with pysam.AlignmentFile(str(output_bam), "wb", header=bam_in.header):
                pass
        return output_bam

    with tempfile.TemporaryDirectory(dir=str(output_bam.parent)) as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        tmp_paths = [str(tmp_dir / f"contig_{index:05d}.bam") for index in range(len(contigs))]
        with ProcessPoolExecutor(max_workers=min(workers, len(contigs))) as pool:
            futures = [
                pool.submit(
                    _reshape_contig_worker,
                    str(source_bam),
                    str(reference_fasta),
                    contig,
                    acceptance_probabilities,
                    seed,
                    min_insert_size,
                    max_insert_size,
                    min_mapq,
                    tmp_path,
                )
                for contig, tmp_path in zip(contigs, tmp_paths)
            ]
            for future in futures:
                future.result()
        pysam.cat("-o", str(output_bam), *tmp_paths)
    return output_bam


@dataclass(frozen=True)
class MixtureReshapeCounts:
    """Per-mixture pair counts from one `count_mixture_reshape_pairs`/
    `reshape_mixture_bam` pass -- the ticket's per-mixture "ā_realised and k"
    bookkeeping, in read-pair units (one count per read1)."""

    tumor_original: int
    tumor_kept: int
    background_original: int
    background_kept: int = 0

    @property
    def realised_acceptance(self) -> float:
        """This mixture's own realised ā -- may differ slightly from the
        source-wide `mean_acceptance` since it depends on the particular
        (downsampled) fragment-length composition of this one mixture."""
        return self.tumor_kept / self.tumor_original if self.tumor_original else 0.0

    @property
    def total_kept(self) -> int:
        return self.tumor_kept + self.background_kept

    def to_dict(self) -> dict[str, object]:
        return {
            "tumor_original": self.tumor_original,
            "tumor_kept": self.tumor_kept,
            "background_original": self.background_original,
            "background_kept": self.background_kept,
            "realised_acceptance": self.realised_acceptance,
            "total_kept": self.total_kept,
        }


def count_mixture_reshape_pairs(
    bam_path: Path,
    reference_fasta: Path,
    tumor_rg_ids: set[str],
    acceptance_table: dict[int, float],
    seed: int,
    min_insert_size: int = DEFAULT_MIN_INSERT_SIZE,
    max_insert_size: int = DEFAULT_MAX_INSERT_SIZE,
    min_mapq: int = DEFAULT_MIN_MAPQ,
) -> MixtureReshapeCounts:
    """First pass over a merged mixture BAM: counts tumor-RG pairs kept under
    `acceptance_table` and background-RG pairs (untouched at this stage), one
    count per read1 -- without writing an output BAM. Used to size the common
    background-thinning target `k` across a batch of mixtures before the
    (expensive) write pass runs `reshape_mixture_bam` on each.

    `min_mapq` (see `_meets_mapq_threshold`) is applied only to tumor-RG
    pairs, matching `reshape_mixture_bam`: the mixture's real cfDNA
    background reads are not a rejection-sampling artefact source (the
    coverage-bias report's finding is specific to tlen-only thinning of the
    mechanically-sheared HCC1395 tumor source) and are left untouched here,
    as `background_original` already documents.
    """
    tumor_original = tumor_kept = background_original = 0
    with pysam.AlignmentFile(str(bam_path), "rb", reference_filename=str(reference_fasta)) as bam:
        for read in bam.fetch(until_eof=True):
            if not _passes_pair_filter(read):
                continue
            if not (read.flag & READ1_FLAG):
                continue
            if read.get_tag("RG") in tumor_rg_ids:
                tumor_original += 1
                if not _meets_mapq_threshold(read, min_mapq):
                    continue
                tlen = abs(read.template_length)
                if tlen == 0 or not (min_insert_size <= tlen <= max_insert_size):
                    continue
                probability = acceptance_table.get(tlen, 0.0)
                if probability <= 0:
                    continue
                if probability >= 1.0 or _acceptance_draw(seed, read.query_name, tlen) < probability:
                    tumor_kept += 1
            else:
                background_original += 1
    return MixtureReshapeCounts(tumor_original, tumor_kept, background_original)


def compute_background_keep_threshold(
    bam_path: Path,
    reference_fasta: Path,
    tumor_rg_ids: set[str],
    background_keep_count: int,
    seed: int,
) -> float:
    """The `background_thin_draw` cutoff `theta` such that keeping every
    background read1 with `draw <= theta` retains *exactly*
    `background_keep_count` background pairs (clamped to
    `[0, background_original]`) -- an exact top-k selection by draw value
    rather than an independent per-pair probability threshold, so the batch's
    common total-pair target (`reshape_admixture_mixture.py`'s
    `common_target_total`) is hit exactly rather than only in expectation.
    Draws are floats in `[0, 1)` derived from a cryptographic hash, so ties
    are astronomically unlikely.

    `background_keep_count <= 0` returns a threshold below every possible
    draw (nothing kept); `background_keep_count` at or above the background
    population returns a threshold above every possible draw (everything
    kept) without needing to know that population size up front.
    """
    if background_keep_count <= 0:
        return -1.0
    draws: list[float] = []
    with pysam.AlignmentFile(str(bam_path), "rb", reference_filename=str(reference_fasta)) as bam:
        for read in bam.fetch(until_eof=True):
            if not _passes_pair_filter(read):
                continue
            if not (read.flag & READ1_FLAG):
                continue
            if read.get_tag("RG") in tumor_rg_ids:
                continue
            draws.append(background_thin_draw(seed, read.query_name))
    if background_keep_count >= len(draws):
        return 1.0
    draws.sort()
    return draws[background_keep_count - 1]


def reshape_mixture_bam(
    bam_path: Path,
    reference_fasta: Path,
    output_bam: Path,
    tumor_rg_ids: set[str],
    acceptance_table: dict[int, float],
    background_keep_threshold: float,
    seed: int,
    min_insert_size: int = DEFAULT_MIN_INSERT_SIZE,
    max_insert_size: int = DEFAULT_MAX_INSERT_SIZE,
    min_mapq: int = DEFAULT_MIN_MAPQ,
) -> tuple[Path, MixtureReshapeCounts]:
    """Writes a reshaped mixture BAM: tumor-RG pairs are thinned by
    `acceptance_table` exactly as `reshape_bam` does (same `_acceptance_draw`
    keyed on `(seed, qname, tlen)`) and additionally dropped below `min_mapq`
    (see `_meets_mapq_threshold`); background-RG pairs are kept iff
    `background_thin_draw(seed, qname) <= background_keep_threshold` (from
    `compute_background_keep_threshold`, an exact top-k cutoff rather than a
    probability), so both mates always agree and the batch's common
    total-pair target is hit exactly. Non-tumor, non-background RGs (there
    should be none in an `admixture.py` mixture) are dropped, matching
    `reshape_bam`'s "reshaping can suppress but never invent" convention.

    `min_mapq` is applied to tumor-RG reads only -- the mixture's background
    reads are left untouched, same reasoning as `count_mixture_reshape_pairs`.

    A single pass, mirroring `reshape_bam`: the input is already
    coordinate-sorted and a subset of a sorted sequence stays sorted, so no
    re-sort is needed before indexing.
    """
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    tumor_original = tumor_kept = background_original = background_kept = 0
    with pysam.AlignmentFile(str(bam_path), "rb", reference_filename=str(reference_fasta)) as bam_in:
        with pysam.AlignmentFile(str(output_bam), "wb", header=bam_in.header) as bam_out:
            for read in bam_in.fetch(until_eof=True):
                if not _passes_pair_filter(read):
                    continue
                is_read1 = bool(read.flag & READ1_FLAG)
                if read.get_tag("RG") in tumor_rg_ids:
                    if is_read1:
                        tumor_original += 1
                    if not _meets_mapq_threshold(read, min_mapq):
                        continue
                    tlen = abs(read.template_length)
                    if tlen == 0 or not (min_insert_size <= tlen <= max_insert_size):
                        continue
                    probability = acceptance_table.get(tlen, 0.0)
                    if probability <= 0:
                        continue
                    if probability < 1.0 and _acceptance_draw(seed, read.query_name, tlen) >= probability:
                        continue
                    if is_read1:
                        tumor_kept += 1
                    bam_out.write(read)
                else:
                    if is_read1:
                        background_original += 1
                    if background_thin_draw(seed, read.query_name) > background_keep_threshold:
                        continue
                    if is_read1:
                        background_kept += 1
                    bam_out.write(read)
    counts = MixtureReshapeCounts(tumor_original, tumor_kept, background_original, background_kept)
    return output_bam, counts


def ks_statistic_from_histograms(histogram_a: Counter, histogram_b: Counter) -> float:
    """Two-sample Kolmogorov-Smirnov statistic computed directly from
    histograms -- the max absolute gap between the two empirical CDFs.

    Equivalent to `scipy.stats.ks_2samp`'s statistic on the underlying
    per-read samples, but never materialises them: a real BAM's fragment
    count (10^8-10^9) would dwarf the distinct-insert-size histogram (10^3),
    same memory argument as `fragment_extraction.compute_fragment_stats`.
    """
    total_a = sum(histogram_a.values())
    total_b = sum(histogram_b.values())
    if total_a == 0 or total_b == 0:
        raise ValueError("cannot compute a KS statistic against an empty histogram")

    support = sorted(set(histogram_a) | set(histogram_b))
    cumulative_a = cumulative_b = 0.0
    max_gap = 0.0
    for insert_size in support:
        cumulative_a += histogram_a.get(insert_size, 0) / total_a
        cumulative_b += histogram_b.get(insert_size, 0) / total_b
        max_gap = max(max_gap, abs(cumulative_a - cumulative_b))
    return max_gap


def _histogram_median(histogram: Counter) -> float:
    """`statistics.median` semantics computed over a histogram -- mirrors
    `fragment_extraction.FragmentExtractor._histogram_median`."""
    total = sum(histogram.values())
    if total == 0:
        raise ValueError("empty histogram has no median")
    lower_rank, upper_rank = (total - 1) // 2, total // 2
    cumulative = 0
    lower_value = None
    for insert_size in sorted(histogram):
        cumulative += histogram[insert_size]
        if lower_value is None and cumulative > lower_rank:
            lower_value = insert_size
        if cumulative > upper_rank:
            return (lower_value + insert_size) / 2
    raise ValueError("empty histogram has no median")


def _histogram_mode(histogram: Counter) -> int:
    if not histogram:
        raise ValueError("empty histogram has no mode")
    peak = max(histogram.values())
    return min(insert_size for insert_size, count in histogram.items() if count == peak)


@dataclass(frozen=True)
class ReshapingVerificationResult:
    statistic: float
    threshold: float
    passed: bool
    reshaped_median: float
    reshaped_mode: int
    target_median: float
    target_mode: int
    n_reshaped: int
    n_target: int

    def to_dict(self) -> dict[str, object]:
        return {
            "statistic": self.statistic,
            "threshold": self.threshold,
            "passed": self.passed,
            "reshaped_median": self.reshaped_median,
            "reshaped_mode": self.reshaped_mode,
            "target_median": self.target_median,
            "target_mode": self.target_mode,
            "n_reshaped": self.n_reshaped,
            "n_target": self.n_target,
        }


def verify_reshaping(
    reshaped_histogram: Counter, target_histogram: Counter, threshold: float
) -> ReshapingVerificationResult:
    statistic = ks_statistic_from_histograms(reshaped_histogram, target_histogram)
    return ReshapingVerificationResult(
        statistic=statistic,
        threshold=threshold,
        passed=statistic <= threshold,
        reshaped_median=_histogram_median(reshaped_histogram),
        reshaped_mode=_histogram_mode(reshaped_histogram),
        target_median=_histogram_median(target_histogram),
        target_mode=_histogram_mode(target_histogram),
        n_reshaped=sum(reshaped_histogram.values()),
        n_target=sum(target_histogram.values()),
    )


@dataclass(frozen=True)
class ReshapingConfig:
    source_bam: Path
    reference_fasta: Path
    output_dir: Path
    # Exactly one of these must be set: `target_bam` computes the target
    # histogram fresh from a background BAM/CRAM's own observed fragment
    # sizes (the ticket's "don't hardcode a shape"); `target_histogram_path`
    # reuses a histogram a previous run already wrote via `write_histogram_csv`
    # (e.g. once per background library, shared across many reshaping runs).
    target_bam: Path | None = None
    target_histogram_path: Path | None = None
    min_insert_size: int = DEFAULT_MIN_INSERT_SIZE
    max_insert_size: int = DEFAULT_MAX_INSERT_SIZE
    ks_threshold: float = DEFAULT_KS_THRESHOLD
    seed: int = 42
    sample_name: str = "reshaped"
    min_mapq: int = DEFAULT_MIN_MAPQ
    # Per-contig process-pool fan-out for the histogram and retention passes
    # (see `compute_insert_size_histogram`/`reshape_bam`). Default 1
    # preserves the original single-threaded, single-pass behavior exactly.
    workers: int = 1

    def __post_init__(self) -> None:
        if (self.target_bam is None) == (self.target_histogram_path is None):
            raise ValueError(
                "ReshapingConfig requires exactly one of target_bam or target_histogram_path, got "
                f"target_bam={self.target_bam!r}, target_histogram_path={self.target_histogram_path!r}"
            )

    @classmethod
    def from_sources(
        cls,
        source_bam: Path,
        reference_fasta: Path,
        output_dir: Path,
        target_bam: Path | None = None,
        target_histogram_path: Path | None = None,
        min_insert_size: int = DEFAULT_MIN_INSERT_SIZE,
        max_insert_size: int = DEFAULT_MAX_INSERT_SIZE,
        ks_threshold: float = DEFAULT_KS_THRESHOLD,
        seed: int = 42,
        sample_name: str = "reshaped",
        min_mapq: int = DEFAULT_MIN_MAPQ,
        workers: int = 1,
    ) -> "ReshapingConfig":
        return cls(
            source_bam=Path(source_bam),
            reference_fasta=Path(reference_fasta),
            output_dir=Path(output_dir),
            target_bam=Path(target_bam) if target_bam is not None else None,
            target_histogram_path=Path(target_histogram_path) if target_histogram_path is not None else None,
            min_insert_size=int(min_insert_size),
            max_insert_size=int(max_insert_size),
            ks_threshold=float(ks_threshold),
            seed=int(seed),
            sample_name=sample_name,
            min_mapq=int(min_mapq),
            workers=int(workers),
        )

    @property
    def output_bam(self) -> Path:
        return self.output_dir / f"{self.sample_name}.reshaped.bam"

    @property
    def verification_report_path(self) -> Path:
        return self.output_dir / f"{self.sample_name}.reshaping_verification.json"


class FragmentReshaper:
    """Drives histogram -> acceptance-probability -> retention -> verification
    for one tumor source BAM, per the ticket's post-alignment design."""

    def __init__(self, config: ReshapingConfig) -> None:
        self.config = config

    def load_target_histogram(self) -> Counter:
        cfg = self.config
        if cfg.target_histogram_path is not None:
            return read_histogram_csv(cfg.target_histogram_path)
        return compute_insert_size_histogram(
            cfg.target_bam, cfg.reference_fasta, cfg.min_insert_size, cfg.max_insert_size, cfg.min_mapq, cfg.workers
        )

    def compute_source_histogram(self) -> Counter:
        cfg = self.config
        return compute_insert_size_histogram(
            cfg.source_bam, cfg.reference_fasta, cfg.min_insert_size, cfg.max_insert_size, cfg.min_mapq, cfg.workers
        )

    def compute_acceptance_probabilities(self, source_histogram: Counter, target_histogram: Counter) -> dict[int, float]:
        return acceptance_probabilities(source_histogram, target_histogram)

    def reshape(self, probabilities: dict[int, float]) -> Path:
        cfg = self.config
        cfg.output_dir.mkdir(parents=True, exist_ok=True)
        reshape_bam(
            cfg.source_bam,
            cfg.reference_fasta,
            cfg.output_bam,
            probabilities,
            cfg.seed,
            cfg.min_insert_size,
            cfg.max_insert_size,
            cfg.min_mapq,
            cfg.workers,
        )
        pysam.index(str(cfg.output_bam))
        return cfg.output_bam

    def verify(self, reshaped_histogram: Counter, target_histogram: Counter) -> ReshapingVerificationResult:
        return verify_reshaping(reshaped_histogram, target_histogram, self.config.ks_threshold)

    def write_verification_report(self, result: ReshapingVerificationResult) -> Path:
        cfg = self.config
        cfg.output_dir.mkdir(parents=True, exist_ok=True)
        cfg.verification_report_path.write_text(json.dumps(result.to_dict(), indent=2))
        return cfg.verification_report_path

    def run(self, raise_on_failure: bool = True) -> dict[str, object]:
        """Runs all three stages and writes the reshaped BAM, its index, the
        three histograms (source/target/reshaped, for later inspection), and
        the verification report.

        Raises when verification fails and `raise_on_failure` (the default) --
        matching `AdmixtureBuilder.verify_downsample_depth`'s convention of
        failing loud rather than silently handing a bad output downstream.
        """
        cfg = self.config
        cfg.output_dir.mkdir(parents=True, exist_ok=True)

        target_histogram = self.load_target_histogram()
        source_histogram = self.compute_source_histogram()
        probabilities = self.compute_acceptance_probabilities(source_histogram, target_histogram)

        output_bam = self.reshape(probabilities)
        reshaped_histogram = compute_insert_size_histogram(
            output_bam, cfg.reference_fasta, cfg.min_insert_size, cfg.max_insert_size, cfg.min_mapq, cfg.workers
        )
        result = self.verify(reshaped_histogram, target_histogram)
        report_path = self.write_verification_report(result)

        write_histogram_csv(source_histogram, cfg.output_dir / f"{cfg.sample_name}.source_histogram.csv")
        write_histogram_csv(target_histogram, cfg.output_dir / f"{cfg.sample_name}.target_histogram.csv")
        write_histogram_csv(reshaped_histogram, cfg.output_dir / f"{cfg.sample_name}.reshaped_histogram.csv")

        if raise_on_failure and not result.passed:
            raise ValueError(
                f"reshaped insert-size histogram for {cfg.sample_name!r} does not match the target "
                f"within threshold: KS statistic {result.statistic:.4f} > {result.threshold:.4f} "
                f"(see {report_path})"
            )

        return {
            "output_bam": output_bam,
            "verification": result,
            "verification_report_path": report_path,
            "n_source_fragments": sum(source_histogram.values()),
            "n_reshaped_fragments": sum(reshaped_histogram.values()),
        }
