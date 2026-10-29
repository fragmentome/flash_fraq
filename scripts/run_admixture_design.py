"""Driver for the admixture feasibility table + mixture plan (phase2 §7 step 1).

Reads the native-coverage manifest and writes admixture_feasibility.csv,
admixture_plan.csv and split_assignment.csv under --output-root. Every knob that
shapes the design -- the f ladder, the candidate target coverages, the retention
floors, the split fractions, the seed -- comes from
workflow/config/admixture_design.yaml via the admixture_design rule, not from
defaults here, so the plan is reproducible from the config alone.
"""

from __future__ import annotations

import json
from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.admixture_design import AdmixtureDesignConfig, AdmixtureDesigner


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--coverage-manifest-path", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--summary-path", type=Path)
    parser.add_argument("--tumor-sample-id", required=True)
    parser.add_argument("--tumor-subject-id")
    parser.add_argument("--tumor-platform", required=True)
    parser.add_argument("--tumor-native-coverage", required=True, type=float)
    parser.add_argument("--tumor-bam-path", required=True)
    parser.add_argument("--matched-pair-sample-id", required=True)
    parser.add_argument("--matched-pair-subject-id")
    parser.add_argument("--matched-pair-source-group", default="hcc1395bl")
    parser.add_argument("--matched-pair-platform", required=True)
    parser.add_argument("--matched-pair-native-coverage", required=True, type=float)
    parser.add_argument("--matched-pair-bam-path", required=True)
    parser.add_argument("--matched-pair-replicates-per-fraction", type=int, default=1)
    parser.add_argument("--fraction-ladder", required=True, type=float, nargs="+")
    parser.add_argument("--candidate-target-coverages", required=True, type=float, nargs="+")
    parser.add_argument("--background-source-group", action="append", default=[])
    parser.add_argument("--background-tumor-status")
    parser.add_argument("--split-fraction", action="append", required=True, metavar="NAME=WEIGHT")
    parser.add_argument("--mixtures-per-fraction", type=int, required=True)
    parser.add_argument("--min-background-samples-per-source", type=int, required=True)
    parser.add_argument("--min-background-retention-per-source", type=float, required=True)
    parser.add_argument("--reference-build", default="hg19")
    parser.add_argument("--seed", type=int, required=True)
    args = parser.parse_args()

    designer = AdmixtureDesigner(
        AdmixtureDesignConfig.from_sources(
            coverage_manifest_path=args.coverage_manifest_path,
            output_root=args.output_root,
            tumor={
                "sample_id": args.tumor_sample_id,
                "subject_id": args.tumor_subject_id or args.tumor_sample_id,
                "platform": args.tumor_platform,
                "native_coverage": args.tumor_native_coverage,
                "bam_path": args.tumor_bam_path,
            },
            matched_pair_background={
                "sample_id": args.matched_pair_sample_id,
                "subject_id": args.matched_pair_subject_id or args.matched_pair_sample_id,
                "source_group": args.matched_pair_source_group,
                "platform": args.matched_pair_platform,
                "native_coverage": args.matched_pair_native_coverage,
                "bam_path": args.matched_pair_bam_path,
            },
            fraction_ladder=args.fraction_ladder,
            candidate_target_coverages=args.candidate_target_coverages,
            split_fractions=dict(_parse_split_fraction(item) for item in args.split_fraction),
            background_source_groups=args.background_source_group,
            background_tumor_status=args.background_tumor_status,
            mixtures_per_fraction=args.mixtures_per_fraction,
            matched_pair_replicates_per_fraction=args.matched_pair_replicates_per_fraction,
            min_background_samples_per_source=args.min_background_samples_per_source,
            min_background_retention_per_source=args.min_background_retention_per_source,
            reference_build=args.reference_build,
            seed=args.seed,
        )
    )
    summary = designer.run()
    if args.summary_path is not None:
        designer.write_summary(summary, args.summary_path)

    print(f"admixture plan: {designer.config.plan_path}")
    print(json.dumps(summary, indent=2))


def _parse_split_fraction(argument: str) -> tuple[str, float]:
    name, separator, weight = argument.partition("=")
    if not separator or not name:
        raise ValueError(f"--split-fraction expects NAME=WEIGHT, got {argument!r}")
    return name, float(weight)


if __name__ == "__main__":
    main()
