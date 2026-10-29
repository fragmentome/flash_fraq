"""Driver for the cohort feature matrix over the real Snyder + Cristiano samples.

Pivots and joins every completed per-sample local_features output under
`--derived-root` into one table, tagged by source dataset/platform. This is the
pre-admixture stand-in for the phase2 §9 matrix -- see flash_fraq.cohort_matrix
for what it deliberately does not carry.

`--dataset NAME=PLATFORM` (repeatable) both selects the datasets to include and
supplies the platform label that the §6 decoy classifier tries to predict.
Invoked by the cohort_feature_matrix rule in workflow/rules/cohort_qc.smk.
"""

from __future__ import annotations

import json
from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.cohort_matrix import CohortFeatureMatrixAssembler, CohortMatrixConfig


def _parse_dataset(argument: str) -> tuple[str, str]:
    name, separator, platform = argument.partition("=")
    if not separator or not name or not platform:
        raise ValueError(f"--dataset expects NAME=PLATFORM, got {argument!r}")
    return name, platform


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--derived-root", required=True, type=Path)
    parser.add_argument("--dataset", required=True, action="append", metavar="NAME=PLATFORM")
    parser.add_argument("--output-path", required=True, type=Path)
    parser.add_argument("--summary-path", type=Path)
    parser.add_argument("--reference-build", default="hg19")
    args = parser.parse_args()

    assembler = CohortFeatureMatrixAssembler(
        CohortMatrixConfig.from_sources(
            derived_root=args.derived_root,
            dataset_platforms=dict(_parse_dataset(dataset) for dataset in args.dataset),
            output_path=args.output_path,
            reference_build=args.reference_build,
        )
    )
    matrix = assembler.assemble()
    output_path = assembler.write(matrix)
    summary = assembler.summarise(matrix)

    if args.summary_path is not None:
        args.summary_path.parent.mkdir(parents=True, exist_ok=True)
        args.summary_path.write_text(json.dumps(summary, indent=2))

    print(f"cohort feature matrix: {output_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
