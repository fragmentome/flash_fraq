"""Driver for MAT-23: align one HCC1395/HCC1395BL FASTQ pair to hg19 and
verify the resulting BAM's @SQ header (docs/phase2-preprocessing-plan.md §0/§1b).

Invoked per-run by workflow/rules/hcc1395.smk's `align_hcc1395` rule.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.alignment import AlignmentConfig, FastqAligner
from flash_fraq.input_manifest import ManifestEntry
from flash_fraq.notifications import notify_discord


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--run-accession", required=True)
    parser.add_argument("--source-group", required=True, choices=("hcc1395", "hcc1395bl"))
    parser.add_argument("--tumor-status", required=True, choices=("tumor", "healthy"))
    parser.add_argument("--platform", required=True)
    parser.add_argument("--fastq-r1", required=True, type=Path)
    parser.add_argument("--fastq-r2", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    # samtools sort takes -m as a PER-THREAD cap, so the sort's peak footprint
    # is roughly threads x this value: the cluster's 8 x 16G fits its 160G
    # mem_mb request, while the 32-thread workstation profile has to shrink it
    # to stay inside 125G of RAM. Deployment-specific, hence configurable
    # (workflow/config/hcc1395_local.yaml) rather than hardcoded here.
    parser.add_argument("--sort-memory-per-thread", default="16G")
    args = parser.parse_args()

    entry = ManifestEntry(
        sample_id=args.run_accession,
        subject_id=args.run_accession,
        source_group=args.source_group,
        tumor_status=args.tumor_status,
        platform=args.platform,
        file_type="fastq",
        fastq_r1_path=args.fastq_r1,
        fastq_r2_path=args.fastq_r2,
    )

    try:
        aligner = FastqAligner(
            AlignmentConfig.from_sources(
                reference_fasta=args.reference_fasta,
                output_dir=args.output_dir,
                threads=args.threads,
                sort_memory_per_thread=args.sort_memory_per_thread,
            )
        )
        bam_path = aligner.align(entry)
        header_matches = aligner.verify_matches_reference_fasta(bam_path, args.reference_fasta)
    except Exception as e:
        notify_discord(f"MAT-23 HCC1395 alignment FAILED for run {args.run_accession}: {e}")
        raise

    print(f"aligned BAM: {bam_path}")
    print(f"@SQ header matches {args.reference_fasta}: {header_matches}")
    if not header_matches:
        notify_discord(
            f"MAT-23 HCC1395 alignment FAILED for run {args.run_accession}: @SQ header does not match {args.reference_fasta}"
        )
        raise SystemExit(1)

    # Alignment is verified above -- the raw FASTQ (freshly downloaded from
    # SRA, not archival like the Cristiano subset) is re-fetchable and no
    # longer needed once the BAM it produced is confirmed good.
    args.fastq_r1.unlink()
    args.fastq_r2.unlink()

    notify_discord(f"MAT-23 HCC1395 alignment finished for run {args.run_accession}: {bam_path}")


if __name__ == "__main__":
    main()
