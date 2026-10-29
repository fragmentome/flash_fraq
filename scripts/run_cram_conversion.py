"""Driver for MAT-24: convert one hg19-aligned BAM (MAT-19 Cristiano or MAT-20
Snyder) to a verified, indexed CRAM, then delete the raw BAM and its .bai.

Pure lossless format conversion -- no deduplication (docs/pipeline-modules.md:28;
dedup stays inside admixture.py's merge_and_dedup(), once, post-merge).

Runs on Institute HPC (isilon-attached) -- invoked per-sample by the
cram_cristiano rule in workflow/rules/mat19.smk and the cram_snyder rule in
workflow/rules/mat20.smk (see workflow/Snakefile).
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.alignment import AlignmentConfig, FastqAligner
from flash_fraq.notifications import notify_discord


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--bam-path", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    # Threads for the samtools view -C encode and the two verification
    # flagstat passes; the rules pass their own {threads} so the value matches
    # what SLURM allocated as cpus-per-task.
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    aligner = FastqAligner(
        AlignmentConfig.from_sources(
            reference_fasta=args.reference_fasta,
            output_dir=args.bam_path.parent,
            threads=args.threads,
        )
    )
    try:
        cram_path = aligner.convert_bam_to_verified_cram(args.bam_path, args.reference_fasta)
    except Exception as e:
        notify_discord(f"MAT-24 CRAM conversion FAILED for {args.bam_path}: {e}")
        raise

    print(f"CRAM: {cram_path}")

    # Conversion already succeeded and was verified above -- a Discord hiccup
    # here must not discard that result.
    try:
        notify_discord(f"MAT-24 CRAM conversion finished: {cram_path}")
    except Exception as e:
        print(f"WARNING: failed to send Discord notification: {e}")


if __name__ == "__main__":
    main()
