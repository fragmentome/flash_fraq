"""Driver for the HCC1395 lane -> library CRAM merge (MAT-23 Block C).

Merges the banked per-run CRAMs of one SEQC2 library
(`<prep>_<site>_<arm>_<variant>`) into a single indexed, quickchecked
sample-level CRAM. See `flash_fraq.hcc1395_merge` for why the library -- not
the cell line and not the site/arm -- is the merge unit.

Invoked by `rule merge_hcc1395_library` in workflow/rules/hcc1395.smk;
`--list-units` is the same planner the rule uses to build its DAG, handy for
checking what would be merged before anything runs.
"""

from __future__ import annotations

import json
from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.hcc1395_merge import (
    DEFAULT_MERGE_PREPS,
    HCC1395_DATASET,
    CramLibraryMerger,
    HCC1395MergeConfig,
    LibraryMergePlanner,
)


def _build_parser() -> ArgumentParser:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--pairs-manifest", required=True, type=Path)
    parser.add_argument("--cram-root", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--dataset", default=HCC1395_DATASET)
    parser.add_argument("--prep", action="append", dest="preps", metavar="PREP")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--verify-read-counts",
        action="store_true",
        help="also flagstat every lane and the merge and require the totals to agree "
        "(a second full pass over hundreds of GB -- off by default)",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--sample-id", help="merged library to build, e.g. HCC1395_WGS_IL_2")
    group.add_argument(
        "--list-units",
        action="store_true",
        help="print the planned merge units as JSON and exit without merging",
    )
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    config = HCC1395MergeConfig.from_sources(
        pairs_manifest=args.pairs_manifest,
        cram_root=args.cram_root,
        reference_fasta=args.reference_fasta,
        dataset=args.dataset,
        preps=tuple(args.preps) if args.preps else DEFAULT_MERGE_PREPS,
        threads=args.threads,
        verify_read_counts=args.verify_read_counts,
    )
    planner = LibraryMergePlanner(config)
    units = planner.plan_available()

    if args.list_units:
        print(
            json.dumps(
                [
                    {
                        "sample_id": unit.sample_id,
                        "library_name": unit.library_name,
                        "subject_id": unit.subject_id,
                        "source_group": unit.source_group,
                        "tumor_status": unit.tumor_status,
                        "runs": list(unit.run_accessions),
                        "approximate_coverage": round(unit.approximate_coverage(), 2),
                        "merged_cram": str(unit.merged_cram_path(config.cram_root, config.dataset)),
                    }
                    for unit in units
                ],
                indent=2,
            )
        )
        return

    unit = next((candidate for candidate in units if candidate.sample_id == args.sample_id), None)
    if unit is None:
        raise SystemExit(
            f"{args.sample_id!r} is not a mergeable library under {args.cram_root}/{args.dataset}; "
            f"available: {[candidate.sample_id for candidate in units]}"
        )

    merged_cram = CramLibraryMerger(config).merge(unit)
    print(
        json.dumps(
            {
                "sample_id": unit.sample_id,
                "library_name": unit.library_name,
                "n_lanes": len(unit.lanes),
                "runs": list(unit.run_accessions),
                "merged_cram": str(merged_cram),
                "approximate_coverage": round(unit.approximate_coverage(), 2),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
