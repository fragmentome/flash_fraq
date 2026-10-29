"""Input manifest schema for cfDNA source files feeding the pipeline (docs/phase2-preprocessing-plan.md §1).

Distinguishes already-aligned BAM sources (§1a) from raw-FASTQ sources needing
alignment (§1b) so downstream stages (admixture.py, cnv_binning.py) can select the
right entry path per sample. Note: the Cristiano DELFI pancreatic subset is
FASTQ-only on the current archive mount (`HumGenTempData`), contradicting
phase1/phase2's assumption of an already-aligned BAM subset at the `ARCHIVE` mount —
see the open Linear ticket tracking reconciliation of this against the docs.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from pandas import DataFrame, read_csv

SOURCE_GROUPS: tuple[str, ...] = ("snyder", "cristiano_delfi", "hcc1395", "hcc1395bl", "1000g_nygc")
TUMOR_STATUSES: tuple[str, ...] = ("tumor", "healthy", "negative_control")
FILE_TYPES: tuple[str, ...] = ("bam", "fastq")


@dataclass(frozen=True)
class ManifestEntry:
    sample_id: str
    subject_id: str
    source_group: str
    tumor_status: str
    platform: str
    file_type: str
    reference_build: str = "hg19"
    reference_build_verified: bool = False
    native_coverage: float | None = None
    bam_path: Path | None = None
    fastq_r1_path: Path | None = None
    fastq_r2_path: Path | None = None

    def __post_init__(self) -> None:
        if self.file_type == "bam" and self.bam_path is None:
            raise ValueError(f"{self.sample_id}: file_type='bam' requires bam_path")
        if self.file_type == "fastq" and (self.fastq_r1_path is None or self.fastq_r2_path is None):
            raise ValueError(f"{self.sample_id}: file_type='fastq' requires fastq_r1_path and fastq_r2_path")


@dataclass(frozen=True)
class InputManifestConfig:
    output_path: Path
    cristiano_fastq_root: Path | None = None

    @classmethod
    def from_sources(cls, output_path: Path, cristiano_fastq_root: Path | None = None) -> "InputManifestConfig":
        return cls(output_path=output_path, cristiano_fastq_root=cristiano_fastq_root)


class InputManifestBuilder:
    """Collects ManifestEntry rows per source group and writes/reads the manifest CSV."""

    def __init__(self, config: InputManifestConfig) -> None:
        self.config = config
        self.entries: list[ManifestEntry] = []

    def discover_cristiano_fastq(self) -> list[ManifestEntry]:
        root = self.config.cristiano_fastq_root
        if root is None:
            raise ValueError("cristiano_fastq_root not set on InputManifestConfig")

        discovered = []
        for cancer_type_dir in sorted(p for p in root.glob("*") if p.is_dir()):
            for r1_path in sorted(cancer_type_dir.glob("*_R1.fastq.gz")):
                sample_id = r1_path.name.removesuffix("_R1.fastq.gz")
                r2_path = cancer_type_dir / f"{sample_id}_R2.fastq.gz"
                if not r2_path.exists():
                    continue
                discovered.append(
                    ManifestEntry(
                        sample_id=sample_id,
                        subject_id=sample_id,
                        source_group="cristiano_delfi",
                        tumor_status="tumor",
                        platform="HiSeq2500",
                        file_type="fastq",
                        reference_build="hg19",
                        reference_build_verified=False,
                        fastq_r1_path=r1_path,
                        fastq_r2_path=r2_path,
                    )
                )
        self.entries.extend(discovered)
        return discovered

    def add_entry(self, entry: ManifestEntry) -> None:
        self.entries.append(entry)

    def write(self) -> Path:
        output_path = self.config.output_path
        output_path.parent.mkdir(parents=True, exist_ok=True)
        DataFrame(
            [
                {
                    "sample_id": e.sample_id,
                    "subject_id": e.subject_id,
                    "source_group": e.source_group,
                    "tumor_status": e.tumor_status,
                    "platform": e.platform,
                    "file_type": e.file_type,
                    "reference_build": e.reference_build,
                    "reference_build_verified": e.reference_build_verified,
                    "native_coverage": e.native_coverage,
                    "bam_path": str(e.bam_path) if e.bam_path is not None else "",
                    "fastq_r1_path": str(e.fastq_r1_path) if e.fastq_r1_path is not None else "",
                    "fastq_r2_path": str(e.fastq_r2_path) if e.fastq_r2_path is not None else "",
                }
                for e in self.entries
            ]
        ).to_csv(output_path, index=False)
        return output_path

    @staticmethod
    def load(manifest_path: Path) -> list[ManifestEntry]:
        rows = read_csv(manifest_path, keep_default_na=False)
        return [
            ManifestEntry(
                sample_id=row["sample_id"],
                subject_id=row["subject_id"],
                source_group=row["source_group"],
                tumor_status=row["tumor_status"],
                platform=row["platform"],
                file_type=row["file_type"],
                reference_build=row["reference_build"],
                reference_build_verified=bool(row["reference_build_verified"]),
                native_coverage=float(row["native_coverage"]) if row["native_coverage"] != "" else None,
                bam_path=Path(row["bam_path"]) if row["bam_path"] else None,
                fastq_r1_path=Path(row["fastq_r1_path"]) if row["fastq_r1_path"] else None,
                fastq_r2_path=Path(row["fastq_r2_path"]) if row["fastq_r2_path"] else None,
            )
            for _, row in rows.iterrows()
        ]
