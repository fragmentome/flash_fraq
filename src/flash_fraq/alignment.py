"""FASTQ alignment against hg19 and post-alignment header verification
(docs/phase2-preprocessing-plan.md §0, §1b; MAT-19, MAT-20).

Used as the raw-FASTQ entry path for sources without a ready aligned BAM
(Cristiano DELFI pancreatic subset, freshly-downloaded Snyder samples).
Header verification compares against an already-aligned hg19 BAM from the
same source family rather than a hardcoded contig list, since the exact
reference FASTA variant (decoy set, contig naming) behind the existing
Snyder/Cristiano BAMs is still an open item (phase2-preprocessing-plan.md §0)
— a build-name match alone isn't sufficient, contig sets must match exactly
for samtools merge to work downstream.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from subprocess import PIPE, CalledProcessError, Popen, run

from flash_fraq.input_manifest import ManifestEntry


@dataclass(frozen=True)
class AlignmentConfig:
    reference_fasta: Path
    output_dir: Path
    threads: int = 8
    sort_memory_per_thread: str = "16G"
    staging_dir: Path | None = None

    @classmethod
    def from_sources(
        cls,
        reference_fasta: Path,
        output_dir: Path,
        threads: int = 8,
        sort_memory_per_thread: str = "16G",
        staging_dir: Path | None = None,
    ) -> "AlignmentConfig":
        return cls(
            reference_fasta=reference_fasta,
            output_dir=output_dir,
            threads=threads,
            sort_memory_per_thread=sort_memory_per_thread,
            staging_dir=staging_dir,
        )


class FastqAligner:
    """BWA-MEM alignment of a FASTQ pair against hg19, plus @SQ header
    verification against a known-good hg19 BAM (phase2-preprocessing-plan.md §0/§1b)."""

    def __init__(self, config: AlignmentConfig) -> None:
        self.config = config

    def align(self, entry: ManifestEntry) -> Path:
        if entry.file_type != "fastq":
            raise ValueError(f"{entry.sample_id}: align() requires file_type='fastq', got {entry.file_type!r}")

        cfg = self.config
        # staging_dir routes bwa/samtools I/O through fast local scratch instead
        # of isilon -- multiple concurrent jobs streaming multi-GB BAMs over NFS
        # otherwise collapse per-job throughput (MAT-20 samples were taking
        # 30-44h of "real time" in bwa mem alone). Only the final indexed BAM
        # is copied out to isilon, below.
        work_dir = cfg.staging_dir if cfg.staging_dir is not None else cfg.output_dir
        work_dir.mkdir(parents=True, exist_ok=True)
        sorted_bam = work_dir / f"{entry.sample_id}.sorted.bam"
        read_group = f"@RG\\tID:{entry.sample_id}\\tSM:{entry.sample_id}\\tPL:{entry.platform}"

        unsorted_bam = work_dir / f"{entry.sample_id}.unsorted.bam"

        bwa_process = Popen(
            [
                "bwa", "mem",
                "-t", str(cfg.threads),
                "-R", read_group,
                str(cfg.reference_fasta),
                str(entry.fastq_r1_path),
                str(entry.fastq_r2_path),
            ],
            stdout=PIPE,
        )
        # samtools 0.1.19 on the cluster doesn't auto-detect SAM on stdin the
        # way modern samtools sort does -- it expects BAM binary, so bwa
        # mem's SAM text has to be converted via `samtools view -bS` first.
        with unsorted_bam.open("wb") as unsorted_fh:
            view_process = Popen(
                ["samtools", "view", "-bS", "-"],
                stdin=bwa_process.stdout,
                stdout=unsorted_fh,
            )
            bwa_process.stdout.close()
            view_process.wait()
        bwa_process.wait()
        if bwa_process.returncode != 0:
            raise CalledProcessError(bwa_process.returncode, "bwa mem")
        if view_process.returncode != 0:
            raise CalledProcessError(view_process.returncode, "samtools view")

        # samtools 1.24 (unlike the old 0.1.19) rejects the classic
        # positional `sort <in.bam> <out.prefix>` form outright ("Use -T
        # PREFIX / -o FILE") -- explicit -o is required.
        run(
            ["samtools", "sort", "-@", str(cfg.threads), "-m", cfg.sort_memory_per_thread, "-o", str(sorted_bam), str(unsorted_bam)],
            check=True,
        )
        unsorted_bam.unlink()

        if not sorted_bam.exists() or sorted_bam.stat().st_size == 0:
            raise RuntimeError(f"samtools sort did not produce {sorted_bam}")

        run(["samtools", "index", str(sorted_bam)], check=True)

        if cfg.staging_dir is None:
            return sorted_bam

        cfg.output_dir.mkdir(parents=True, exist_ok=True)
        final_bam = cfg.output_dir / f"{entry.sample_id}.sorted.bam"
        shutil.copy2(sorted_bam, final_bam)
        shutil.copy2(f"{sorted_bam}.bai", f"{final_bam}.bai")
        shutil.rmtree(cfg.staging_dir)
        return final_bam

    def get_header_contigs(self, bam_path: Path, reference_fasta: Path | None = None) -> dict[str, int]:
        cmd = ["samtools", "view", "-H"]
        if reference_fasta is not None:
            cmd += ["--reference", str(reference_fasta)]
        cmd.append(str(bam_path))
        result = run(cmd, check=True, capture_output=True, text=True)
        contigs = {}
        for line in result.stdout.splitlines():
            if not line.startswith("@SQ"):
                continue
            tags = dict(tag.split(":", 1) for tag in line.split("\t")[1:])
            contigs[tags["SN"]] = int(tags["LN"])
        return contigs

    def verify_matches_reference_header(self, bam_path: Path, reference_bam_path: Path) -> bool:
        return self.get_header_contigs(bam_path) == self.get_header_contigs(reference_bam_path)

    def get_fasta_contigs(self, reference_fasta_path: Path) -> dict[str, int]:
        fai_path = Path(f"{reference_fasta_path}.fai")
        contigs = {}
        with open(fai_path) as handle:
            for line in handle:
                name, length, *_ = line.split("\t")
                contigs[name] = int(length)
        return contigs

    def verify_matches_reference_fasta(self, bam_path: Path, reference_fasta_path: Path) -> bool:
        return self.get_header_contigs(bam_path) == self.get_fasta_contigs(reference_fasta_path)

    def convert_to_cram(self, bam_path: Path, reference_fasta: Path) -> Path:
        """Losslessly convert an indexed, sorted BAM to an indexed CRAM
        (MAT-24). Pure format conversion -- no deduplication."""
        cram_path = bam_path.with_suffix(".cram")
        # -@: CRAM encoding is the slowest step of the conversion and is the
        # one part of it that threads well -- measured 19 MB/s at one thread
        # against 40 MB/s at eight. On the HCC1395 cohort, whose sorted BAMs
        # run 130-350 GB, that is the difference between ~2 h and ~1 h per
        # median run, and it is time the raw BAM spends occupying a bam_align
        # slot on the fileset.
        run(
            [
                "samtools", "view", "-C", "-@", str(self.config.threads),
                "--reference", str(reference_fasta),
                "-o", str(cram_path), str(bam_path),
            ],
            check=True,
        )
        run(["samtools", "index", "-@", str(self.config.threads), str(cram_path)], check=True)
        return cram_path

    def get_flagstat_total_reads(self, alignment_path: Path, reference_fasta: Path | None = None) -> int:
        # flagstat has no top-level --reference flag (unlike view/merge) --
        # CRAM decoding is instead selected via --input-fmt-option.
        cmd = ["samtools", "flagstat", "-@", str(self.config.threads)]
        if reference_fasta is not None:
            cmd += ["--input-fmt-option", f"reference={reference_fasta}"]
        cmd.append(str(alignment_path))
        result = run(cmd, check=True, capture_output=True, text=True)
        return int(result.stdout.splitlines()[0].split()[0])

    def verify_cram_conversion(self, bam_path: Path, cram_path: Path, reference_fasta: Path) -> bool:
        """Confirm a CRAM produced by convert_to_cram() is a lossless copy of
        bam_path: matching @SQ header and an exact read count match (MAT-24) --
        this is pure format conversion, so an exact match is required, not
        just a non-zero read count."""
        headers_match = self.get_header_contigs(bam_path) == self.get_header_contigs(cram_path, reference_fasta)
        bam_reads = self.get_flagstat_total_reads(bam_path)
        cram_reads = self.get_flagstat_total_reads(cram_path, reference_fasta)
        return headers_match and bam_reads == cram_reads

    def convert_bam_to_verified_cram(self, bam_path: Path, reference_fasta: Path) -> Path:
        """Convert bam_path to an indexed CRAM, verify it losslessly
        reproduces the BAM, then delete the raw BAM and its .bai (MAT-24).
        Raises RuntimeError -- and leaves the BAM in place -- if verification
        fails."""
        cram_path = self.convert_to_cram(bam_path, reference_fasta)
        if not self.verify_cram_conversion(bam_path, cram_path, reference_fasta):
            raise RuntimeError(f"CRAM verification failed for {cram_path} against {bam_path}")

        bam_path.unlink()
        bai_path = Path(f"{bam_path}.bai")
        if bai_path.exists():
            bai_path.unlink()
        return cram_path
