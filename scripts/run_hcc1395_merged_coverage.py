"""Collect the merged HCC1395 libraries' mosdepth summaries into one CSV
(MAT-23 Block C, task 3).

`rule hcc1395_merged_coverage` runs mosdepth with the exact flags
`rule localfeat_coverage` uses, so these numbers sit on the same scale as the
547 rows already in `native_coverage.csv`. This script reads the `mean` column
of each summary's `total` row -- via
`NativeCoverageHarvester.parse_mosdepth_summary`, the same parser the cohort
harvest uses -- and emits the table plus the exact lines to paste over the
`native_coverage: 35.0` / `34.0` placeholders in
`workflow/config/admixture_design.yaml`.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.hcc1395_merge import (
    DEFAULT_MERGE_PREPS,
    HCC1395_DATASET,
    HCC1395MergeConfig,
    LibraryMergePlanner,
    MergedCoverageReporter,
)


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--pairs-manifest", required=True, type=Path)
    parser.add_argument("--cram-root", required=True, type=Path)
    parser.add_argument("--derived-root", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--dataset", default=HCC1395_DATASET)
    parser.add_argument("--prep", action="append", dest="preps", metavar="PREP")
    parser.add_argument("--output-path", required=True, type=Path)
    parser.add_argument(
        "--pairs-only",
        action="store_true",
        default=True,
        help="Report only libraries in a complete tumor/normal pair (default; "
        "matches workflow/config/hcc1395_merge.yaml's pairs_only and the set "
        "rule hcc1395_merge_library actually merges).",
    )
    parser.add_argument("--no-pairs-only", action="store_false", dest="pairs_only")
    args = parser.parse_args()

    config = HCC1395MergeConfig.from_sources(
        pairs_manifest=args.pairs_manifest,
        cram_root=args.cram_root,
        reference_fasta=args.reference_fasta,
        dataset=args.dataset,
        preps=tuple(args.preps) if args.preps else DEFAULT_MERGE_PREPS,
    )
    planner = LibraryMergePlanner(config)
    units = planner.plan_paired_available() if args.pairs_only else planner.plan_available()
    if not units:
        raise SystemExit(f"no mergeable HCC1395 libraries under {args.cram_root}/{args.dataset}")

    reporter = MergedCoverageReporter(config)
    rows = reporter.rows(units, args.derived_root)
    output_path = reporter.write(rows, args.output_path)

    print(f"merged-library native coverage: {output_path}")
    print(reporter.config_yaml_block(rows))


if __name__ == "__main__":
    main()
