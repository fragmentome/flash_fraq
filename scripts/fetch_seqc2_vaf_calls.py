"""Download the SEQC2 high-confidence somatic SNV VAF call set for
HCC1395/HCC1395BL (docs/phase1-sourcing.md), needed by MAT-13's admixture
empirical-VAF fraction check.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path
from subprocess import run

SEQC2_VAF_CALLS_URL: str = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/"


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    run(
        ["wget", "-r", "-np", "-nH", "--cut-dirs=5", "-P", str(args.output_dir), SEQC2_VAF_CALLS_URL],
        check=True,
    )
    print(f"downloaded SEQC2 somatic VAF call set to {args.output_dir}")


if __name__ == "__main__":
    main()
