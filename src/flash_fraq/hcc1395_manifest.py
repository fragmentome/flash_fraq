"""Resolve per-sample SRR accessions for HCC1395/HCC1395BL (SEQC2) from SRA
study SRP162370 via the SRA Run Selector runinfo endpoint
(docs/phase1-sourcing.md; MAT-23).

Filters the study-wide runinfo table down to Illumina WGS runs for
HCC1395/HCC1395BL specifically, rather than a blind `prefetch SRP162370`
full-project pull -- the study also contains Ion Torrent, PacBio, and
Nanopore runs plus derived/replicate samples outside the tumor/normal WGS
scope this project needs.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

import requests

from flash_fraq.input_manifest import ManifestEntry

SRA_RUNINFO_URL: str = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"
HCC1395_SRA_STUDY_ACCESSION: str = "SRP162370"


@dataclass(frozen=True)
class SraRun:
    run_accession: str
    sample_name: str
    center_name: str
    platform: str
    library_strategy: str

    @property
    def source_group(self) -> str:
        return "hcc1395bl" if "BL" in self.sample_name.upper() else "hcc1395"

    @property
    def tumor_status(self) -> str:
        return "healthy" if self.source_group == "hcc1395bl" else "tumor"


@dataclass(frozen=True)
class HCC1395ManifestConfig:
    output_path: Path
    sra_study_accession: str = HCC1395_SRA_STUDY_ACCESSION

    @classmethod
    def from_sources(
        cls, output_path: Path, sra_study_accession: str = HCC1395_SRA_STUDY_ACCESSION
    ) -> "HCC1395ManifestConfig":
        return cls(output_path=output_path, sra_study_accession=sra_study_accession)


class HCC1395RunResolver:
    """Resolves SRR run accessions for HCC1395/HCC1395BL WGS runs from the
    SRA Run Selector runinfo table (docs/phase1-sourcing.md)."""

    def __init__(self, config: HCC1395ManifestConfig) -> None:
        self.config = config

    def fetch_runinfo_csv(self) -> str:
        response = requests.get(
            SRA_RUNINFO_URL,
            params={"acc": self.config.sra_study_accession},
            timeout=60,
        )
        response.raise_for_status()
        return response.text

    @staticmethod
    def parse_runinfo(runinfo_csv: str) -> list[SraRun]:
        reader = csv.DictReader(StringIO(runinfo_csv))
        return [
            SraRun(
                run_accession=row["Run"],
                sample_name=row["SampleName"],
                center_name=row["CenterName"],
                platform=row["Platform"],
                library_strategy=row["LibraryStrategy"],
            )
            for row in reader
        ]

    def resolve_wgs_illumina_runs(self) -> list[SraRun]:
        runs = self.parse_runinfo(self.fetch_runinfo_csv())
        return [
            run
            for run in runs
            if run.library_strategy == "WGS"
            and run.platform.upper() == "ILLUMINA"
            and run.sample_name.upper().startswith("HCC1395")
        ]

    def to_manifest_entries(self, runs: list[SraRun], fastq_root: Path) -> list[ManifestEntry]:
        entries = []
        for run in runs:
            sample_dir = fastq_root / run.run_accession
            entries.append(
                ManifestEntry(
                    sample_id=run.run_accession,
                    subject_id=run.sample_name,
                    source_group=run.source_group,
                    tumor_status=run.tumor_status,
                    platform=run.platform,
                    file_type="fastq",
                    fastq_r1_path=sample_dir / f"{run.run_accession}_1.fastq.gz",
                    fastq_r2_path=sample_dir / f"{run.run_accession}_2.fastq.gz",
                )
            )
        return entries

    def write_run_list(self) -> Path:
        runs = self.resolve_wgs_illumina_runs()
        output_path = self.config.output_path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["run_accession", "sample_name", "center_name", "platform", "source_group", "tumor_status"])
            for run in runs:
                writer.writerow(
                    [run.run_accession, run.sample_name, run.center_name, run.platform, run.source_group, run.tumor_status]
                )
        return output_path
