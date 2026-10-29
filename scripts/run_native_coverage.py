"""Driver for the native-coverage harvest (phase2 §7 step 1 inputs).

Walks the local_features derived tree, reads the `total` row of every completed
sample's mosdepth summary, and writes it as an input_manifest CSV with
`native_coverage` populated -- the number `AdmixtureConfig` requires and that
nothing else in the repo produces.

`--dataset` is repeatable and takes NAME:PLATFORM:SOURCE_GROUP:TUMOR_STATUS; it
both selects the dataset directories to walk and supplies the cohort metadata a
mosdepth summary cannot carry. `--phenotype NAME=PATH` gives a dataset a per-sample
`tumor_status` table (manifests/cristiano_delfi_phenotype.csv for Cristiano). Invoked by the native_coverage_manifest rule in
workflow/rules/admixture_design.smk.
"""

from __future__ import annotations

import json
from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.admixture_design import NativeCoverageConfig, NativeCoverageHarvester

DATASET_FIELDS = ("name", "platform", "source_group", "tumor_status")


def _parse_dataset(argument: str) -> dict[str, str]:
    parts = argument.split(":")
    if len(parts) != len(DATASET_FIELDS) or not all(parts):
        raise ValueError(f"--dataset expects {':'.join(field.upper() for field in DATASET_FIELDS)}, got {argument!r}")
    return dict(zip(DATASET_FIELDS, parts))


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--derived-root", required=True, type=Path)
    parser.add_argument("--cram-root", required=True, type=Path)
    parser.add_argument("--dataset", required=True, action="append", metavar="NAME:PLATFORM:GROUP:STATUS")
    parser.add_argument("--phenotype", action="append", default=[], metavar="NAME=PATH")
    parser.add_argument("--output-path", required=True, type=Path)
    parser.add_argument("--summary-path", type=Path)
    parser.add_argument("--reference-build", default="hg19")
    args = parser.parse_args()
    phenotype_paths = dict(item.split("=", 1) for item in args.phenotype)
    datasets = [_parse_dataset(dataset) for dataset in args.dataset]
    unknown = set(phenotype_paths) - {dataset["name"] for dataset in datasets}
    if unknown:
        raise SystemExit(f"--phenotype names datasets not given with --dataset: {sorted(unknown)}")
    for dataset in datasets:
        dataset["phenotype_path"] = phenotype_paths.get(dataset["name"])

    harvester = NativeCoverageHarvester(
        NativeCoverageConfig.from_sources(
            derived_root=args.derived_root,
            cram_root=args.cram_root,
            datasets=datasets,
            output_path=args.output_path,
            reference_build=args.reference_build,
        )
    )
    entries = harvester.harvest()
    if not entries:
        raise SystemExit(f"no completed samples with a mosdepth summary under {args.derived_root}")
    output_path = harvester.write(entries)
    summary = harvester.summarise(entries)

    if args.summary_path is not None:
        args.summary_path.parent.mkdir(parents=True, exist_ok=True)
        args.summary_path.write_text(json.dumps(summary, indent=2))

    print(f"native coverage manifest: {output_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
