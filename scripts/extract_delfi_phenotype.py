"""Per-sample phenotype for the Cristiano DELFI release, read off the CRAM headers.

The EGA release ships no phenotype table with the EGAF accessions, but the
alignment did: every CRAM's `bwa mem` @PG line names the FASTQ it was aligned
from, and those FASTQs were downloaded into one directory per DELFI cohort
(`.../FASTQ/healthy/EGAF..._R1.fastq.gz`, `.../FASTQ/lung_cancer/...`). That
directory name is the phenotype. `cohort` keeps it verbatim; `tumor_status` maps
it onto input_manifest's vocabulary (`healthy` / `tumor`).

Writes manifests/cristiano_delfi_phenotype.csv. Fails if any sample's header
does not name exactly one cohort directory, so the table never has gaps.
"""

from __future__ import annotations

import csv
import re
import subprocess
from argparse import ArgumentParser
from pathlib import Path

HEALTHY_COHORT = "healthy"
PHENOTYPE_COLUMNS = ("sample_id", "cohort", "tumor_status", "source")


def cohort_from_header(header: str, sample_id: str) -> str:
    bwa_lines = [line for line in header.splitlines() if line.startswith("@PG") and "\tID:bwa" in line]
    cohorts = {match for line in bwa_lines for match in re.findall(rf"/FASTQ/([^/\s]+)/{sample_id}_R1", line)}
    if len(cohorts) != 1:
        raise ValueError(f"{sample_id}: expected one FASTQ cohort directory in the bwa @PG line, found {sorted(cohorts)}")
    return cohorts.pop()


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cram-dir", type=Path, default=Path("/mnt/8tb_hdd/flash_fraq_crams/Cristiano_Velculescu_DELFI/CRAM")
    )
    parser.add_argument("--output-path", type=Path, default=Path("manifests/cristiano_delfi_phenotype.csv"))
    args = parser.parse_args()

    sample_ids = (args.cram_dir / "sample_list.txt").read_text().split()
    rows = []
    for sample_id in sorted(sample_ids):
        cram = args.cram_dir / sample_id / f"{sample_id}.sorted.cram"
        header = subprocess.run(["samtools", "view", "-H", str(cram)], check=True, capture_output=True, text=True).stdout
        cohort = cohort_from_header(header, sample_id)
        rows.append(
            {
                "sample_id": sample_id,
                "cohort": cohort,
                "tumor_status": "healthy" if cohort == HEALTHY_COHORT else "tumor",
                "source": "cram_header_bwa_fastq_dir",
            }
        )

    args.output_path.parent.mkdir(parents=True, exist_ok=True)
    with args.output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PHENOTYPE_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} samples -> {args.output_path}")


if __name__ == "__main__":
    main()
