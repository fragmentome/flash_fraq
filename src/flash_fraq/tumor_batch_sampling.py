"""One-pass multi-output tumor downsample sampler (docs/phase2-preprocessing-plan.md §7 perf fix).

`AdmixtureBuilder.subsample_bam`'s `samtools view -s SEED.FRACTION` has no
index-assisted subsampling path -- every invocation does a full linear decode
of the source CRAM/BAM regardless of how small a fraction it keeps. The
admixture plan's rows (both the `real_background` and `matched_pair` regimes
-- see `admixture_design.py`'s single global `TumorSource`) all downsample
from the *same* tumor CRAM, so building each row's tumor_sub.bam independently
means decoding that one ~88GB file once per row (~490 times for the current
plan).

`TumorBatchSampler` decodes it once and writes every row's tumor_sub.bam in
that same pass. It cannot reuse samtools' own per-read hash (there is no
public API for it), so it substitutes a different deterministic per-read draw
(`zlib.crc32` on the query name -- stable across runs and processes, unlike
Python's salted built-in `hash()`) combined with a per-target phase offset
derived from each row's seed (a Weyl/low-discrepancy sequence:
`(seed * golden_ratio_conjugate) % 1.0`), so that one shared draw per read
yields an independent-looking inclusion decision per target fraction. Both
mates of a pair share `query_name` and therefore get the same decision.

This does NOT reproduce samtools -s's exact selection -- a row built by this
sampler will not contain bit-identical reads to the same row built by the old
per-row `samtools view -s` path. What MAT-13 actually needs is preserved:
a fixed seed always retains the same reads (rerun stability), and distinct
seeds retain different, independent-looking subsets (no two mixtures
degenerately sharing their whole retained set). Already-completed mixture
outputs built under the old path are untouched and are not expected to match
newly-built ones bit-for-bit.
"""

from __future__ import annotations

import zlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pysam

_GOLDEN_RATIO_CONJUGATE = 0.6180339887498949
DEFAULT_BATCH_SIZE = 50_000


@dataclass(frozen=True)
class TumorDownsampleTarget:
    """One admixture_plan.csv row's tumor-side downsample output."""

    output_bam: Path
    seed: int
    fraction: float

    def __post_init__(self) -> None:
        if not 0 < self.fraction <= 1:
            raise ValueError(f"fraction must be in (0, 1], got {self.fraction}")

    @property
    def phase(self) -> float:
        return (self.seed * _GOLDEN_RATIO_CONJUGATE) % 1.0


def _read_draws(query_names: list[str]) -> np.ndarray:
    return np.array([zlib.crc32(name.encode()) / 2**32 for name in query_names], dtype=np.float64)


class TumorBatchSampler:
    """Decodes a tumor BAM/CRAM exactly once, writing every target's downsampled BAM."""

    def __init__(self, batch_size: int = DEFAULT_BATCH_SIZE) -> None:
        self.batch_size = batch_size

    def run(
        self,
        tumor_bam: Path,
        reference_fasta: Path,
        targets: list[TumorDownsampleTarget],
        threads: int = 1,
    ) -> int:
        """Write every `target.output_bam`; return the total read count in `tumor_bam`.

        The total is a free byproduct of the single pass this method already
        makes. Callers should use it in place of a per-row `samtools flagstat`
        on the full source (the `reads_before` half of
        `AdmixtureBuilder.verify_downsample_depth`), which would otherwise
        reintroduce the same ~490x redundant full-file decode this sampler
        exists to eliminate -- flagstat still has to walk every record.
        """
        if not targets:
            return 0

        phases = np.array([target.phase for target in targets], dtype=np.float64)
        fractions = np.array([target.fraction for target in targets], dtype=np.float64)

        for target in targets:
            target.output_bam.parent.mkdir(parents=True, exist_ok=True)

        with pysam.AlignmentFile(
            str(tumor_bam), "r", reference_filename=str(reference_fasta), threads=threads
        ) as source:
            writers = [pysam.AlignmentFile(str(target.output_bam), "wb", template=source) for target in targets]
            try:
                total_reads = 0
                batch: list = []
                for read in source:
                    batch.append(read)
                    if len(batch) >= self.batch_size:
                        total_reads += self._flush(batch, phases, fractions, writers)
                        batch = []
                if batch:
                    total_reads += self._flush(batch, phases, fractions, writers)
            finally:
                for writer in writers:
                    writer.close()
        return total_reads

    def _flush(self, batch: list, phases: np.ndarray, fractions: np.ndarray, writers: list) -> int:
        draws = _read_draws([read.query_name for read in batch])
        # (batch_size, n_targets) boolean decision matrix via broadcasting -- the
        # numpy-vectorized replacement for an O(reads x targets) Python loop,
        # which at this source's scale (billions of reads x hundreds of targets)
        # would take on the order of a day rather than the minutes this takes.
        decisions = ((draws[:, None] + phases[None, :]) % 1.0) < fractions[None, :]
        for target_index, writer in enumerate(writers):
            for read_index in np.nonzero(decisions[:, target_index])[0]:
                writer.write(batch[read_index])
        return len(batch)
