"""Driver for the real phase2 §9 feature matrix (MAT-17) over the 495 admixture mixtures.

`flash_fraq.feature_matrix.FeatureMatrixAssembler` expects one aggregated table per
upstream source (all samples' rows already concatenated) under a single input
directory -- that is how MAT-13/14/15/16 write output for a *cohort* run via
Snakemake. The admixture mixtures were instead produced by the three
`scripts/drivers/admixture_mixture_*_driver.sh` drivers, one process per mixture,
each writing into that mixture's own flat `mix_<name>/` directory (see those
drivers' docstrings for why: the mixtures are BAMs in a layout the cohort rules'
CRAM-mirror discovery doesn't cover).

This script bridges the two: it concatenates each of the 495 mixtures'
`sample_metadata.csv`, `fragments/fragment_stats.csv`, `fragments/arm_zscores.csv`,
`fragments/window_ratios.csv`, `ichor/cnv_bins.csv`, and `ichor/ichor_estimates.csv`
into one aggregated table per source, reuses the cohort-wide
`split_assignment.csv` MAT-13's `MixturePlanner` already wrote (subject-level, before
mixing -- shared across the cohort and mixture matrices by construction), writes the
aggregated inputs to `--agg-dir`, then runs the unchanged `FeatureMatrixAssembler`
against them.

`--admixture-plan-path` supplies `flash_fraq.tf_model`'s required `regime` column
(`matched_pair`/`real_background`), joined in on `sample_name` == `sample_id` after
assembly -- `FeatureMatrixAssembler` itself has no notion of regime, only
`admixture_design.MixturePlanner`'s plan does.
"""

from __future__ import annotations

import json
from argparse import ArgumentParser
from pathlib import Path

from pandas import concat, read_csv

from flash_fraq.feature_matrix import FeatureMatrixAssembler, FeatureMatrixConfig
from flash_fraq.genome_mask import DEFAULT_MASK_MARGIN_BP

PER_MIXTURE_SOURCES: dict[str, str] = {
    "sample_metadata.csv": "sample_metadata.csv",
    "fragments/fragment_stats.csv": "fragment_stats.csv",
    "fragments/arm_zscores.csv": "arm_zscores.csv",
    "fragments/window_ratios.csv": "window_ratios.csv",
    "ichor/cnv_bins.csv": "cnv_bins.csv",
    "ichor/ichor_estimates.csv": "ichor_estimates.csv",
}


def _discover_mixture_dirs(mixtures_root: Path) -> list[Path]:
    return sorted(p for p in mixtures_root.iterdir() if p.is_dir() and p.name.startswith("mix_"))


def _aggregate(mixture_dirs: list[Path], relative_path: str) -> tuple[list, list]:
    frames = []
    missing = []
    for mixture_dir in mixture_dirs:
        source_path = mixture_dir / relative_path
        if not source_path.exists():
            missing.append(mixture_dir.name)
            continue
        frames.append(read_csv(source_path))
    return frames, missing


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--mixtures-root", required=True, type=Path)
    parser.add_argument("--split-assignment-path", required=True, type=Path)
    parser.add_argument("--admixture-plan-path", required=True, type=Path)
    parser.add_argument("--agg-dir", required=True, type=Path)
    parser.add_argument("--output-path", required=True, type=Path)
    parser.add_argument("--summary-path", type=Path)
    parser.add_argument(
        "--mask-margin-bp",
        type=int,
        default=DEFAULT_MASK_MARGIN_BP,
        help="Drop cnv_logr_bin_* bins within this many bp of a centromere or "
        "chromosome end before pivoting (0 disables masking).",
    )
    args = parser.parse_args()

    mixture_dirs = _discover_mixture_dirs(args.mixtures_root)
    args.agg_dir.mkdir(parents=True, exist_ok=True)

    missing_by_source: dict[str, list[str]] = {}
    for relative_path, output_name in PER_MIXTURE_SOURCES.items():
        frames, missing = _aggregate(mixture_dirs, relative_path)
        if missing:
            missing_by_source[relative_path] = missing
        concat(frames, ignore_index=True).to_csv(args.agg_dir / output_name, index=False)

    if missing_by_source:
        raise SystemExit(
            f"aborting: missing per-mixture source files: {json.dumps(missing_by_source, indent=2)}"
        )

    split_assignment = read_csv(args.split_assignment_path)
    split_assignment.to_csv(args.agg_dir / "split_assignment.csv", index=False)

    config = FeatureMatrixConfig.from_sources(
        input_dir=args.agg_dir, output_path=args.output_path, mask_margin_bp=args.mask_margin_bp
    )
    assembler = FeatureMatrixAssembler(config)
    matrix = assembler.assemble()

    plan = read_csv(args.admixture_plan_path, usecols=["sample_name", "regime"])
    plan = plan.rename(columns={"sample_name": "sample_id"})
    n_before = matrix.shape[0]
    matrix = matrix.merge(plan, on="sample_id", how="left")
    if matrix["regime"].isna().any():
        unmatched = matrix.loc[matrix["regime"].isna(), "sample_id"].tolist()
        raise SystemExit(f"aborting: {len(unmatched)} sample_id(s) missing from admixture plan: {unmatched}")
    assert matrix.shape[0] == n_before, "regime join must not change row count"

    output_path = assembler.write(matrix)

    summary = {
        "n_mixtures_found": len(mixture_dirs),
        "n_rows": int(matrix.shape[0]),
        "n_columns": int(matrix.shape[1]),
        "primary_feature_columns": len(assembler.get_primary_feature_columns(matrix)),
        "split_group_counts": matrix["split_group"].value_counts().to_dict(),
    }
    if args.summary_path is not None:
        args.summary_path.parent.mkdir(parents=True, exist_ok=True)
        args.summary_path.write_text(json.dumps(summary, indent=2))

    print(f"feature matrix: {output_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
