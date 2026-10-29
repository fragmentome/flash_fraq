"""Driver for one in silico fragment-size reshaping run
(docs/todo/completed/fragment-size-reshaping.md).

Thin CLI over `ReshapingConfig.from_sources` + `FragmentReshaper.run()` --
reads a tumor source BAM/CRAM (HCC1395), retains read pairs by an
insert-size-dependent acceptance probability so the retained population's
fragment-size distribution statistically matches a target ctDNA-like
distribution, and writes the reshaped BAM plus a verification report.

Produces an upstream input that a later wiring step would hand to
`scripts/run_admixture.py` as an alternative `--tumor-bam` for the
`real_background` regime -- that wiring is out of scope here.
"""

from __future__ import annotations

import json
from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.fragment_reshaping import (
    DEFAULT_KS_THRESHOLD,
    DEFAULT_MAX_INSERT_SIZE,
    DEFAULT_MIN_INSERT_SIZE,
    DEFAULT_MIN_MAPQ,
    FragmentReshaper,
    ReshapingConfig,
)


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--source-bam", required=True, type=Path, help="Tumor BAM/CRAM to reshape (e.g. HCC1395).")
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument(
        "--target-bam", type=Path, help="Background BAM/CRAM whose own observed fragment sizes are the reshaping target."
    )
    target_group.add_argument(
        "--target-histogram-path", type=Path, help="Precomputed target histogram CSV (from a previous run's *.target_histogram.csv)."
    )
    parser.add_argument("--min-insert-size", type=int, default=DEFAULT_MIN_INSERT_SIZE)
    parser.add_argument("--max-insert-size", type=int, default=DEFAULT_MAX_INSERT_SIZE)
    parser.add_argument("--ks-threshold", type=float, default=DEFAULT_KS_THRESHOLD)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sample-name", default="reshaped")
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=DEFAULT_MIN_MAPQ,
        help="Minimum MAPQ (both mates) for a tumor pair to be eligible for retention; also applied "
        "when computing the source/target insert-size histograms. See coverage_bias_report.json.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Per-contig process-pool fan-out for the histogram and retention passes (see "
        "fragment_reshaping.compute_insert_size_histogram/reshape_bam). Default 1 preserves the "
        "original single-threaded, single-pass behavior exactly; >1 parallelizes both full-BAM "
        "passes across contigs, writing per-contig BAMs that are then concatenated in header "
        "contig order. Pair accept/reject decisions are deterministic given (seed, qname, insert "
        "size), so the retained set is identical regardless of worker count.",
    )
    args = parser.parse_args()

    reshaper = FragmentReshaper(
        ReshapingConfig.from_sources(
            source_bam=args.source_bam,
            reference_fasta=args.reference_fasta,
            output_dir=args.output_dir,
            target_bam=args.target_bam,
            target_histogram_path=args.target_histogram_path,
            min_insert_size=args.min_insert_size,
            max_insert_size=args.max_insert_size,
            ks_threshold=args.ks_threshold,
            seed=args.seed,
            sample_name=args.sample_name,
            min_mapq=args.min_mapq,
            workers=args.workers,
        )
    )
    result = reshaper.run()
    print(
        json.dumps(
            {
                "sample_name": args.sample_name,
                "output_bam": str(result["output_bam"]),
                "n_source_fragments": result["n_source_fragments"],
                "n_reshaped_fragments": result["n_reshaped_fragments"],
                "verification": result["verification"].to_dict(),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
