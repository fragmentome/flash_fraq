"""Recompute empirical_tumor_fraction for one already-built admixture mixture, in place.

Fixes the fields[9]-only bug in `AdmixtureBuilder.compute_empirical_tumor_fraction`
(see git history / docs/todo for the writeup): the merged/dedup BAM keeps the tumor
and healthy sources' distinct @RG SM tags, so bcftools mpileup emitted one sample
column per source and the old code read only the first one. No re-downsampling or
re-merging is needed -- the dedup BAM on disk is unaffected by that bug, only the
label computed from it was wrong. This script re-runs just
`compute_empirical_tumor_fraction` against the existing dedup BAM for both SNV-VAF
labels, and rewrites `sample_metadata.csv` and the per-label `.txt` files in place.

Also picks up finding 7 of docs/results-corruption-audit.md for free: the estimator
lives in `compute_empirical_tumor_fraction` itself (the depth-weighted ratio of sums,
not the old unweighted mean of per-site ratios), so this script's one mpileup rerun
gets the fields[9] fix and the estimator fix at once, and also writes the per-site
`.sites.tsv.gz` dump the estimator now produces as a side effect.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from pandas import read_csv

from flash_fraq.admixture import AdmixtureBuilder, AdmixtureConfig


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--mixture-dir", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--snv-vaf-path", action="append", required=True, help="label=path, repeatable")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--dedup-bam-path",
        type=Path,
        default=None,
        help=(
            "Override the dedup BAM read for mpileup (e.g. an SSD-staged copy), while "
            "still writing all outputs (sample_metadata.csv, per-label .txt, .sites.tsv.gz) "
            "under --mixture-dir as usual. Defaults to <mixture-dir>/<sample_name>.dedup.bam."
        ),
    )
    args = parser.parse_args()

    snv_vaf_paths = {}
    for item in args.snv_vaf_path:
        label, path = item.split("=", 1)
        snv_vaf_paths[label] = Path(path)

    mixture_dir = args.mixture_dir
    sample_name = mixture_dir.name
    dedup_bam = args.dedup_bam_path if args.dedup_bam_path is not None else mixture_dir / f"{sample_name}.dedup.bam"
    metadata_path = mixture_dir / "sample_metadata.csv"

    metadata = read_csv(metadata_path)
    assert len(metadata) == 1, f"expected exactly one row in {metadata_path}"

    cfg = AdmixtureConfig(
        tumor_bam=Path("/dev/null"),
        healthy_bam=Path("/dev/null"),
        reference_fasta=args.reference_fasta,
        snv_vaf_paths=snv_vaf_paths,
        output_dir=mixture_dir,
        tumor_fraction=0.0,
        target_coverage=1.0,
        tumor_native_coverage=1.0,
        healthy_native_coverage=1.0,
        sample_name=sample_name,
        threads=args.threads,
    )
    builder = AdmixtureBuilder(cfg)

    empirical_tumor_fractions: dict[str, float] = {}
    for label, snv_vaf_path in snv_vaf_paths.items():
        snvs = builder.load_snvs(snv_vaf_path)
        fraction = builder.compute_empirical_tumor_fraction(dedup_bam, snvs, label)
        empirical_tumor_fractions[label] = fraction
        (mixture_dir / f"{sample_name}.empirical_tumor_fraction.{label}.txt").write_text(f"{fraction}\n")

    canonical_label = next(iter(snv_vaf_paths))
    metadata.loc[0, "empirical_tumor_fraction"] = empirical_tumor_fractions[canonical_label]
    for label, fraction in empirical_tumor_fractions.items():
        metadata.loc[0, f"empirical_tumor_fraction_{label}"] = fraction
    metadata.to_csv(metadata_path, index=False)

    print(f"{sample_name}: {empirical_tumor_fractions}")


if __name__ == "__main__":
    main()
