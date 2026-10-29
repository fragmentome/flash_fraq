import re
import subprocess
from collections import Counter
from dataclasses import dataclass
from math import log2
from pathlib import Path

import pysam
from pandas import DataFrame


# samtools view -f 0x2 -F 0x400+0x100+0x800: properly paired, not duplicate,
# not secondary, not supplementary.
REQUIRE_FLAGS: int = 0x2
EXCLUDE_FLAGS: int = 0x400 | 0x100 | 0x800
READ1_FLAG: int = 0x40

# DELFI short/long fragment classification (phase0-spec.md §1.2 / Cristiano et al. 2019).
SHORT_FRAGMENT_RANGE: tuple[int, int] = (100, 150)
LONG_FRAGMENT_RANGE: tuple[int, int] = (151, 220)

WINDOW_INDEX_PATTERN: re.Pattern[str] = re.compile(r"(\d+)$")

# The three per-fragment .bed intermediates run_full()/run() write into
# output_dir. All three are pure intermediates: fragment_midpoints.bed and
# fragments_delfi_5mb.bed are consumed only within this module's own
# run_full() (write_fragment_stats and write_window_fragment_counts
# respectively, both called before run_full() returns); fragments_ichor_1mb.bed
# has no reader anywhere in the codebase -- ichorCNA's readCounter reads the
# BAM directly (see cnv_binning.py), not this file. All three are cheaply
# reproducible from the source BAM + the (static) window bed files, so once
# run_full()'s own outputs (window_fragment_counts.tsv, fragment_stats.csv)
# exist on disk, nothing needs these three files again. See
# cleanup_fragment_intermediates() and docs/results-corruption-audit.md's
# fragment-bed cleanup note.
FRAGMENT_BED_INTERMEDIATES: tuple[str, ...] = (
    "fragment_midpoints.bed",
    "fragments_delfi_5mb.bed",
    "fragments_ichor_1mb.bed",
)


def cleanup_fragment_intermediates(output_dir: Path, *, dry_run: bool = False) -> list[Path]:
    """Removes the per-fragment .bed intermediates from a completed fragment-extraction
    output_dir (see FRAGMENT_BED_INTERMEDIATES for why this is safe).

    Callers are responsible for confirming the downstream outputs that actually need
    these files (window_fragment_counts.tsv within this module, and -- one stage further
    down the chain -- GC correction's window_ratios.csv/arm_zscores.csv) already exist and
    are non-empty before calling this; it does not check that itself, since "downstream
    outputs exist" means different things depending on how far the caller's pipeline has
    progressed. Missing files are skipped silently (idempotent under retries/partial
    cleanup). Returns the list of paths removed (or, if dry_run, that would be removed).
    """
    removed = []
    for name in FRAGMENT_BED_INTERMEDIATES:
        path = output_dir / name
        if not path.exists():
            continue
        removed.append(path)
        if not dry_run:
            path.unlink()
    return removed


@dataclass(frozen=True)
class FragmentExtractionConfig:
    input_bam: Path
    output_dir: Path
    delfi_windows_bed: Path
    ichor_bins_bed: Path
    reference_fasta: Path
    sample_id: str
    bedtools_binary: Path

    @classmethod
    def from_sources(
        cls,
        input_bam: Path,
        output_dir: Path,
        delfi_windows_bed: Path,
        ichor_bins_bed: Path,
        reference_fasta: Path,
        sample_id: str | None = None,
        bedtools_binary: Path = Path("bedtools"),
    ) -> "FragmentExtractionConfig":
        return cls(
            input_bam=input_bam,
            output_dir=output_dir,
            delfi_windows_bed=delfi_windows_bed,
            ichor_bins_bed=ichor_bins_bed,
            reference_fasta=reference_fasta,
            sample_id=sample_id if sample_id is not None else Path(input_bam).stem,
            bedtools_binary=bedtools_binary,
        )


class FragmentExtractor:
    def __init__(self, config: FragmentExtractionConfig) -> None:
        self.config = config

    def run(self) -> tuple[Path, Path]:
        midpoints_bed = self.extract_midpoints()
        return self.intersect_windows(midpoints_bed)

    def run_full(self) -> tuple[Path, Path, Path, Path]:
        """Runs extraction, window intersection, GC/count aggregation, and fragment-stat summary.

        Returns (delfi_intersected_bed, ichor_intersected_bed, window_fragment_counts_path,
        fragment_stats_path).
        """
        midpoints_bed = self.extract_midpoints()
        delfi_out, ichor_out = self.intersect_windows(midpoints_bed)
        window_counts_path = self.write_window_fragment_counts(delfi_out)
        fragment_stats_path = self.write_fragment_stats(midpoints_bed)
        return delfi_out, ichor_out, window_counts_path, fragment_stats_path

    def extract_midpoints(self) -> Path:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        midpoints_bed = self.config.output_dir / "fragment_midpoints.bed"
        # reference_filename is what makes this work on CRAM. htslib otherwise
        # falls back to the @SQ UR: tag, which on these files still points at
        # the cluster path the alignment ran from
        # (/home/isilon/.../hg19.fa) -- absent locally, so decoding fails with
        # a misleading "truncated file" several containers in. Harmless for
        # BAM, which ignores the reference entirely.
        with pysam.AlignmentFile(
            str(self.config.input_bam),
            "rb",
            reference_filename=str(self.config.reference_fasta),
        ) as bam, midpoints_bed.open("w") as out:
            for read in bam.fetch(until_eof=True):
                if not self._passes_filters(read):
                    continue
                tlen = abs(read.template_length)
                if tlen == 0:
                    continue
                fragment_start = min(read.reference_start, read.next_reference_start)
                midpoint = fragment_start + tlen // 2
                out.write(f"{read.reference_name}\t{midpoint}\t{midpoint + 1}\t{tlen}\n")
        return midpoints_bed

    def intersect_windows(self, midpoints_bed: Path) -> tuple[Path, Path]:
        delfi_out = self.config.output_dir / "fragments_delfi_5mb.bed"
        ichor_out = self.config.output_dir / "fragments_ichor_1mb.bed"
        self._run_bedtools_intersect(midpoints_bed, self.config.delfi_windows_bed, delfi_out)
        self._run_bedtools_intersect(midpoints_bed, self.config.ichor_bins_bed, ichor_out)
        return delfi_out, ichor_out

    def aggregate_window_counts(self, delfi_intersected: Path) -> DataFrame:
        """Aggregates delfi-intersected fragments into per-window short/long counts + GC content.

        Short/long classification follows the DELFI convention (100-150bp short,
        151-220bp long); fragments outside both ranges are excluded from the counts.
        """
        counts: dict[int, dict[str, object]] = {}
        with delfi_intersected.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                tlen = int(fields[3])
                window_chrom, window_start, window_end, window_name = fields[4], int(fields[5]), int(fields[6]), fields[7]
                window_index = self._parse_window_index(window_name)

                entry = counts.setdefault(
                    window_index,
                    {"chrom": window_chrom, "start": window_start, "end": window_end, "short_count": 0, "long_count": 0},
                )
                if SHORT_FRAGMENT_RANGE[0] <= tlen <= SHORT_FRAGMENT_RANGE[1]:
                    entry["short_count"] += 1
                elif LONG_FRAGMENT_RANGE[0] <= tlen <= LONG_FRAGMENT_RANGE[1]:
                    entry["long_count"] += 1

        rows = []
        with pysam.FastaFile(str(self.config.reference_fasta)) as fasta:
            for window_index, entry in sorted(counts.items()):
                gc_content = self._window_gc_content(fasta, entry["chrom"], entry["start"], entry["end"])
                rows.append(
                    {
                        "sample_id": self.config.sample_id,
                        "window_index": window_index,
                        "chrom": entry["chrom"],
                        "start": entry["start"],
                        "end": entry["end"],
                        "gc_content": gc_content,
                        "short_count": entry["short_count"],
                        "long_count": entry["long_count"],
                    }
                )
        return DataFrame(rows)

    def write_window_fragment_counts(self, delfi_intersected: Path) -> Path:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.config.output_dir / "window_fragment_counts.tsv"
        self.aggregate_window_counts(delfi_intersected).to_csv(output_path, sep="\t", index=False)
        return output_path

    def compute_fragment_stats(self, midpoints_bed: Path) -> tuple[float, float, float]:
        """Genome-wide fragment_size_median/mode/entropy from the extracted fragment TLENs.

        Streamed into a TLEN histogram rather than materialised as a list.
        fragment_midpoints.bed carries one line per fragment -- tens of GB and
        hundreds of millions of rows on a deep sample -- so read_text() plus a
        per-fragment list would need far more RAM than the rule requests. The
        histogram is bounded by the number of *distinct* fragment lengths
        (order 10^3), and median/mode/entropy are all recoverable from it.
        """
        histogram = self._tlen_histogram(midpoints_bed)
        if not histogram:
            raise ValueError(f"no fragments in {midpoints_bed}; cannot summarise fragment sizes")
        fragment_size_median = self._histogram_median(histogram)
        peak = max(histogram.values())
        fragment_size_mode = min(tlen for tlen, count in histogram.items() if count == peak)
        fragment_size_entropy = self._shannon_entropy(histogram)
        return fragment_size_median, fragment_size_mode, fragment_size_entropy

    @staticmethod
    def _tlen_histogram(midpoints_bed: Path) -> Counter:
        histogram: Counter = Counter()
        with midpoints_bed.open() as handle:
            for line in handle:
                if not line.strip():
                    continue
                histogram[int(line.split("\t")[3])] += 1
        return histogram

    @staticmethod
    def _histogram_median(histogram: Counter) -> float:
        """statistics.median semantics (mean of the two central values on an
        even count) computed over a histogram instead of the raw sequence."""
        total = sum(histogram.values())
        lower_rank, upper_rank = (total - 1) // 2, total // 2
        cumulative = 0
        lower_value = None
        for tlen in sorted(histogram):
            cumulative += histogram[tlen]
            if lower_value is None and cumulative > lower_rank:
                lower_value = tlen
            if cumulative > upper_rank:
                return (lower_value + tlen) / 2
        raise ValueError("empty fragment-length histogram")

    def write_fragment_stats(self, midpoints_bed: Path) -> Path:
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.config.output_dir / "fragment_stats.csv"
        fragment_size_median, fragment_size_mode, fragment_size_entropy = self.compute_fragment_stats(midpoints_bed)
        DataFrame(
            [
                {
                    "sample_id": self.config.sample_id,
                    "fragment_size_median": fragment_size_median,
                    "fragment_size_mode": fragment_size_mode,
                    "fragment_size_entropy": fragment_size_entropy,
                }
            ]
        ).to_csv(output_path, index=False)
        return output_path

    def _parse_window_index(self, window_name: str) -> int:
        match = WINDOW_INDEX_PATTERN.search(window_name)
        if match is None:
            raise ValueError(f"could not parse a numeric window index from {window_name!r}")
        return int(match.group(1))

    def _window_gc_content(self, fasta: pysam.FastaFile, chrom: str, start: int, end: int) -> float:
        sequence = fasta.fetch(chrom, start, end).upper()
        gc_bases = sequence.count("G") + sequence.count("C")
        called_bases = len(sequence) - sequence.count("N")
        if called_bases == 0:
            return float("nan")
        return gc_bases / called_bases

    @staticmethod
    def _shannon_entropy(histogram: Counter) -> float:
        total = sum(histogram.values())
        return -sum((count / total) * log2(count / total) for count in histogram.values())

    def _passes_filters(self, read: pysam.AlignedSegment) -> bool:
        if read.flag & EXCLUDE_FLAGS:
            return False
        if not (read.flag & REQUIRE_FLAGS):
            return False
        if not (read.flag & READ1_FLAG):
            return False
        return True

    def _run_bedtools_intersect(self, midpoints_bed: Path, windows_bed: Path, output_bed: Path) -> None:
        with output_bed.open("w") as out:
            subprocess.run(
                [
                    str(self.config.bedtools_binary),
                    "intersect",
                    "-a", str(midpoints_bed),
                    "-b", str(windows_bed),
                    "-wa", "-wb",
                ],
                stdout=out,
                check=True,
            )
