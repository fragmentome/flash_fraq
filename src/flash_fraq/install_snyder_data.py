"""Download GSE71378 (Snyder) FASTQ from ENA (docs/phase2-preprocessing-plan.md §1a; MAT-20).

Resolves fastq_ftp links per SRX accession and downloads the FASTQ pair for
a sample, producing a ManifestEntry ready for FastqAligner. Supersedes the
old manifest-only script, which resolved ENA links but never downloaded
anything.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from subprocess import run

import requests

from flash_fraq.input_manifest import ManifestEntry

ENA_FILEREPORT_URL: str = "https://www.ebi.ac.uk/ena/portal/api/filereport"

SNYDER_SAMPLES: dict[str, str] = {
    "BH01": "SRX1120757", "IA01": "SRX1120758", "IA02": "SRX1120759", "IA03": "SRX1120760",
    "IA04": "SRX1120761", "IA05": "SRX1120762", "IA06": "SRX1120763", "IA07": "SRX1120764",
    "IA08": "SRX1120765", "IC02": "SRX1120766", "IC03": "SRX1120767", "IC04": "SRX1120768",
    "IC05": "SRX1120769", "IC06": "SRX1120770", "IC08": "SRX1120771", "IC09": "SRX1120772",
    "IC10": "SRX1120773", "IC11": "SRX1120774", "IC12": "SRX1120775", "IC13": "SRX1120776",
    "IC14": "SRX1120777", "IC15": "SRX1120778", "IC16": "SRX1120779", "IC17": "SRX1120780",
    "IC18": "SRX1120781", "IC19": "SRX1120782", "IC20": "SRX1120783", "IC21": "SRX1120784",
    "IC22": "SRX1120785", "IC23": "SRX1120786", "IC24": "SRX1120787", "IC25": "SRX1120788",
    "IC26": "SRX1120789", "IC27": "SRX1120790", "IC28": "SRX1120791", "IC29": "SRX1120792",
    "IC30": "SRX1120793", "IC32": "SRX1120794", "IC33": "SRX1120795", "IC34": "SRX1120796",
    "IC35": "SRX1120797", "IC36": "SRX1120798", "IC37": "SRX1120799", "IC38": "SRX1120800",
    "IC39": "SRX1120801", "IC40": "SRX1120802", "IC41": "SRX1120803", "IC42": "SRX1120804",
    "IC43": "SRX1120805", "IC44": "SRX1120806", "IC46": "SRX1120807", "IC47": "SRX1120808",
    "IC48": "SRX1120809", "IC49": "SRX1120810", "IC50": "SRX1120811", "IC51": "SRX1120812",
    "IC52": "SRX1120813", "IH01": "SRX1120814", "IH02": "SRX1120815", "IH03": "SRX1120816",
}


@dataclass(frozen=True)
class EnaRun:
    title: str
    srx: str
    run_accession: str
    submitted_format: str
    fastq_ftp: tuple[str, ...]


@dataclass(frozen=True)
class SnyderInstallConfig:
    output_dir: Path
    samples: dict[str, str] = field(default_factory=lambda: dict(SNYDER_SAMPLES))

    @classmethod
    def from_sources(cls, output_dir: Path) -> "SnyderInstallConfig":
        return cls(output_dir=output_dir)


class SnyderInstaller:
    """Resolves ENA fastq_ftp links per SRX accession and downloads FASTQ pairs
    for GSE71378 (Snyder) samples (docs/phase2-preprocessing-plan.md §1a)."""

    def __init__(self, config: SnyderInstallConfig) -> None:
        self.config = config

    def resolve_ena_runs(self, title: str, srx: str) -> list[EnaRun]:
        response = requests.get(
            ENA_FILEREPORT_URL,
            params={
                "accession": srx,
                "result": "read_run",
                "fields": "run_accession,submitted_ftp,submitted_format,fastq_ftp",
                "format": "json",
            },
            timeout=30,
        )
        response.raise_for_status()
        return [
            EnaRun(
                title=title,
                srx=srx,
                run_accession=row.get("run_accession", ""),
                submitted_format=row.get("submitted_format", ""),
                fastq_ftp=tuple(row["fastq_ftp"].split(";")) if row.get("fastq_ftp") else (),
            )
            for row in response.json()
        ]

    @staticmethod
    def _split_pair(fastq_ftp: tuple[str, ...]) -> tuple[str, str] | None:
        # ENA filereport commonly lists both an interleaved/combined FASTQ and
        # the split _1/_2 mates for one run (e.g. BH01's SRR2129993 has 3
        # entries) -- pick the split mates specifically rather than assuming
        # exactly 2 entries means a pair.
        r1_url = next((u for u in fastq_ftp if u.endswith("_1.fastq.gz")), None)
        r2_url = next((u for u in fastq_ftp if u.endswith("_2.fastq.gz")), None)
        return (r1_url, r2_url) if r1_url and r2_url else None

    def download_fastq_pair(self, ena_run: EnaRun) -> tuple[Path, Path]:
        pair = self._split_pair(ena_run.fastq_ftp)
        if pair is None:
            raise ValueError(
                f"{ena_run.title} ({ena_run.run_accession}): expected split _1/_2 fastq_ftp entries, got {ena_run.fastq_ftp}"
            )

        sample_dir = self.config.output_dir / ena_run.title
        sample_dir.mkdir(parents=True, exist_ok=True)

        r1_url, r2_url = pair
        r1_path = sample_dir / Path(r1_url).name
        r2_path = sample_dir / Path(r2_url).name
        for url, dest in ((r1_url, r1_path), (r2_url, r2_path)):
            # The same samples kept failing identically across repeated
            # full-cohort resubmissions -- not what a one-off transient blip
            # looks like. Up to `jobs: 30` concurrent install_align_snyder
            # jobs hit ENA's FTP at once, which reads as congestion/rate
            # limiting rather than a single hiccup; -nv (instead of -q)
            # surfaces wget's actual error in the rule log instead of
            # nothing, --random-wait desynchronizes the concurrent retries,
            # and the longer/more numerous backoff gives ENA more room to
            # recover.
            run(
                ["wget", "-nv", "--tries=8", "--waitretry=30", "--random-wait", "-O", str(dest), f"https://{url}"],
                check=True,
            )
        return r1_path, r2_path

    def classify_tumor_status(self, title: str) -> str:
        return "tumor" if title.startswith("IC") else "healthy"

    def install_sample(self, title: str) -> ManifestEntry:
        srx = self.config.samples[title]
        ena_runs = self.resolve_ena_runs(title, srx)
        paired_run = next(ena_run for ena_run in ena_runs if self._split_pair(ena_run.fastq_ftp) is not None)
        fastq_r1, fastq_r2 = self.download_fastq_pair(paired_run)
        return ManifestEntry(
            sample_id=title,
            subject_id=title,
            source_group="snyder",
            tumor_status=self.classify_tumor_status(title),
            platform="HiSeq2000",
            file_type="fastq",
            fastq_r1_path=fastq_r1,
            fastq_r2_path=fastq_r2,
        )
