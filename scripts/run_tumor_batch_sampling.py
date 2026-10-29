"""Driver for the one-pass tumor-side downsample batch (phase2 §7 perf fix).

Reads `admixture_plan.csv` and, for every row, takes the tumor-side downsample
target (output path + seed + `tumor_downsample_fraction`) straight from the
plan -- `admixture_design.MixturePlanner` already computed that arithmetic, so
this script does not re-derive it, same convention as `run_admixture.py`.

Every row's `tumor_bam_path` is the same source CRAM (admixture_design's
single global `TumorSource`, shared by both the `matched_pair` and
`real_background` regimes), so this script decodes it exactly once via
`flash_fraq.tumor_batch_sampling.TumorBatchSampler` instead of the one
independent `samtools view -s` full decode each row used to do.

Invoked by `checkpoint admixture_tumor_downsample_batch` in
workflow/rules/admixture.smk, ahead of every `admixture_mixture` job.
"""

from __future__ import annotations

import csv
import resource
from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.tumor_batch_sampling import TumorBatchSampler, TumorDownsampleTarget

TOTAL_READS_FILENAME = "tumor_total_reads.txt"


def _load_targets(plan_path: Path, output_dir: Path) -> tuple[Path, list[TumorDownsampleTarget]]:
    with open(plan_path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{plan_path}: no rows")

    tumor_bam_paths = {row["tumor_bam_path"] for row in rows}
    if len(tumor_bam_paths) != 1:
        raise ValueError(
            f"{plan_path}: expected exactly one tumor_bam_path across all rows -- a single "
            f"shared tumor source is what makes a one-pass batch sampler correct -- got "
            f"{sorted(tumor_bam_paths)}"
        )
    tumor_bam = Path(next(iter(tumor_bam_paths)))

    # tumor_downsample_fraction == 0 is a legitimate row (a pure-background
    # negative control, nominal admixture fraction f=0 -- see
    # AdmixtureBuilder.run's matching special case): there is no tumor read
    # to keep, so no tumor_sub.bam is needed for it at all, and
    # TumorDownsampleTarget's own `0 < fraction <= 1` validation would reject
    # constructing one anyway. Skip these rows entirely rather than trying to
    # build a (meaningless) empty-fraction target for them.
    targets = [
        TumorDownsampleTarget(
            output_bam=output_dir / f"{row['sample_name']}.tumor_sub.bam",
            seed=int(row["seed"]),
            fraction=float(row["tumor_downsample_fraction"]),
        )
        for row in rows
        if float(row["tumor_downsample_fraction"]) > 0
    ]
    return tumor_bam, targets


def _raise_open_file_limit(minimum: int) -> None:
    """One writer per plan row is open simultaneously -- bump the soft
    RLIMIT_NOFILE if the default (often 1024) would not cover the plan's
    row count, capped at whatever the hard limit already allows."""
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    if soft < minimum:
        resource.setrlimit(resource.RLIMIT_NOFILE, (min(minimum, hard), hard))


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    tumor_bam, targets = _load_targets(args.plan, args.output_dir)
    _raise_open_file_limit(len(targets) + 64)

    sampler = TumorBatchSampler()
    total_reads = sampler.run(tumor_bam, args.reference_fasta, targets, threads=args.threads)

    (args.output_dir / TOTAL_READS_FILENAME).write_text(f"{total_reads}\n")


if __name__ == "__main__":
    main()
