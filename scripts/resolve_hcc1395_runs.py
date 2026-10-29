"""Resolve HCC1395/HCC1395BL SRR run accessions from SRA study SRP162370 and
write a run-list CSV consumed by the Snakemake retrieval/alignment rules
(docs/phase1-sourcing.md; MAT-23).

Invoked by workflow/rules/hcc1395.smk's `resolve_hcc1395_runs` checkpoint.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.hcc1395_manifest import HCC1395_SRA_STUDY_ACCESSION, HCC1395ManifestConfig, HCC1395RunResolver


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--sra-study-accession", default=HCC1395_SRA_STUDY_ACCESSION)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    resolver = HCC1395RunResolver(
        HCC1395ManifestConfig.from_sources(output_path=args.output, sra_study_accession=args.sra_study_accession)
    )
    output_path = resolver.write_run_list()
    print(f"wrote run list: {output_path}")


if __name__ == "__main__":
    main()
