"""Driver for one in silico tumor/healthy admixture mixture (phase2 §7, MAT-13).

Thin CLI over `AdmixtureConfig.from_sources` + `AdmixtureBuilder.run()` -- one
invocation builds one `admixture_plan.csv` row's mixture (downsample, merge,
dedup, both SNV-VAF empirical tumor-fraction estimators, sample_metadata.csv).
Invoked by `rule admixture_mixture` in workflow/rules/admixture.smk, which
supplies every argument straight from a plan row plus the SNV VAF paths and
reference FASTA from workflow/config/admixture.yaml.

`--snv-vaf-path LABEL=PATH` is repeatable so MAT-13's "compute both estimators"
requirement (full-set regression sites, clonal-ratio sites) needs no code
change here to add a third: whatever labels are passed are the labels
recorded in sample_metadata.csv's `empirical_tumor_fraction_<label>` columns.
"""

from __future__ import annotations

import json
from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.admixture import AdmixtureBuilder, AdmixtureConfig


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--tumor-bam", required=True, type=Path)
    parser.add_argument("--healthy-bam", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--snv-vaf-path", action="append", required=True, metavar="LABEL=PATH")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--tumor-fraction", required=True, type=float)
    parser.add_argument("--target-coverage", required=True, type=float)
    parser.add_argument("--tumor-native-coverage", required=True, type=float)
    parser.add_argument("--healthy-native-coverage", required=True, type=float)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sample-name", default="mixed")
    parser.add_argument("--tumor-source-id")
    parser.add_argument("--healthy-source-id")
    parser.add_argument("--tumor-subject-id")
    parser.add_argument("--healthy-subject-id")
    parser.add_argument("--platform-tumor", default="unknown")
    parser.add_argument("--platform-healthy", default="unknown")
    parser.add_argument("--reference-build", default="hg19")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument(
        "--tumor-sub-bam",
        type=Path,
        default=None,
        help=(
            "Precomputed tumor-side downsample (from run_tumor_batch_sampling.py's single "
            "shared-CRAM pass) to use in place of this script's own samtools view -s."
        ),
    )
    parser.add_argument(
        "--tumor-total-reads",
        type=int,
        default=None,
        help=(
            "Precomputed total read count for --tumor-bam (also a byproduct of "
            "run_tumor_batch_sampling.py's pass), to skip a redundant samtools flagstat "
            "over the full tumor source in verify_downsample_depth."
        ),
    )
    args = parser.parse_args()

    builder = AdmixtureBuilder(
        AdmixtureConfig.from_sources(
            tumor_bam=args.tumor_bam,
            healthy_bam=args.healthy_bam,
            reference_fasta=args.reference_fasta,
            snv_vaf_paths=dict(_parse_snv_vaf_path(item) for item in args.snv_vaf_path),
            output_dir=args.output_dir,
            tumor_fraction=args.tumor_fraction,
            target_coverage=args.target_coverage,
            tumor_native_coverage=args.tumor_native_coverage,
            healthy_native_coverage=args.healthy_native_coverage,
            seed=args.seed,
            sample_name=args.sample_name,
            tumor_source_id=args.tumor_source_id,
            healthy_source_id=args.healthy_source_id,
            tumor_subject_id=args.tumor_subject_id,
            healthy_subject_id=args.healthy_subject_id,
            platform_tumor=args.platform_tumor,
            platform_healthy=args.platform_healthy,
            reference_build=args.reference_build,
            threads=args.threads,
        )
    )
    empirical_tumor_fractions = builder.run(
        precomputed_tumor_sub=args.tumor_sub_bam,
        precomputed_tumor_total_reads=args.tumor_total_reads,
    )
    print(json.dumps({"sample_name": args.sample_name, **empirical_tumor_fractions}, indent=2))


def _parse_snv_vaf_path(argument: str) -> tuple[str, Path]:
    label, separator, path = argument.partition("=")
    if not separator or not label:
        raise ValueError(f"--snv-vaf-path expects LABEL=PATH, got {argument!r}")
    return label, Path(path)


if __name__ == "__main__":
    main()
