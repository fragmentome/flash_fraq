"""Driver for MAT-19: align one Cristiano DELFI pancreatic FASTQ pair to hg19
and verify the resulting BAM's @SQ header, as a small-subset check before
aligning the full sample set (docs/phase2-preprocessing-plan.md §1b).

Runs on Institute HPC (isilon-attached) — invoked per-sample by the
align_cristiano rule in workflow/rules/mat19.smk (see workflow/Snakefile).
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.alignment import AlignmentConfig, FastqAligner
from flash_fraq.input_manifest import InputManifestConfig, InputManifestBuilder
from flash_fraq.notifications import notify_discord


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True, help="Cristiano sample_id to align (e.g. an EGA run accession)")
    parser.add_argument("--cristiano-fastq-root", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    try:
        manifest_builder = InputManifestBuilder(
            InputManifestConfig.from_sources(
                output_path=args.output_dir / "cristiano_manifest.csv",
                cristiano_fastq_root=args.cristiano_fastq_root,
            )
        )
        entries = manifest_builder.discover_cristiano_fastq()
        entry = next((e for e in entries if e.sample_id == args.sample_id), None)
        if entry is None:
            available = ", ".join(e.sample_id for e in entries)
            raise SystemExit(f"sample_id {args.sample_id!r} not found under {args.cristiano_fastq_root}. Available: {available}")

        aligner = FastqAligner(
            AlignmentConfig.from_sources(
                reference_fasta=args.reference_fasta,
                output_dir=args.output_dir,
                threads=args.threads,
            )
        )
        bam_path = aligner.align(entry)
        header_matches = aligner.verify_matches_reference_fasta(bam_path, args.reference_fasta)
    except Exception as e:
        notify_discord(f"MAT-19 alignment FAILED for sample {args.sample_id}: {e}")
        raise

    print(f"aligned BAM: {bam_path}")
    print(f"@SQ header matches {args.reference_fasta}: {header_matches}")
    if not header_matches:
        notify_discord(f"MAT-19 alignment FAILED for sample {args.sample_id}: @SQ header does not match {args.reference_fasta}")
        raise SystemExit(1)

    notify_discord(f"MAT-19 alignment finished for sample {args.sample_id}: {bam_path}")


if __name__ == "__main__":
    main()
