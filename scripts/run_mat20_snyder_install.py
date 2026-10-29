"""Driver for MAT-20: download and align one GSE71378 (Snyder) sample fresh
from ENA, verify the resulting BAM's @SQ header, as a small-subset check
before extending to the full 60-sample set (docs/phase2-preprocessing-plan.md §1a).

Runs on Institute HPC (isilon-attached) — invoked per-sample by the
install_align_snyder rule in workflow/rules/mat20.smk (see workflow/Snakefile).
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.alignment import AlignmentConfig, FastqAligner
from flash_fraq.install_snyder_data import SnyderInstallConfig, SnyderInstaller
from flash_fraq.notifications import notify_discord


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--sample-title", default="BH01", help="Snyder sample title, e.g. BH01")
    parser.add_argument("--download-dir", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--staging-dir",
        type=Path,
        default=None,
        help="Local scratch root to download and align on instead of --download-dir/--output-dir directly; "
        "only the final BAM+BAI are copied out to --output-dir",
    )
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--sort-memory-per-thread", default="4G")
    args = parser.parse_args()

    snyder_output_dir = args.staging_dir if args.staging_dir is not None else args.download_dir
    alignment_staging_dir = args.staging_dir / args.sample_title if args.staging_dir is not None else None

    try:
        installer = SnyderInstaller(SnyderInstallConfig.from_sources(output_dir=snyder_output_dir))
        entry = installer.install_sample(args.sample_title)

        aligner = FastqAligner(
            AlignmentConfig.from_sources(
                reference_fasta=args.reference_fasta,
                output_dir=args.output_dir,
                threads=args.threads,
                sort_memory_per_thread=args.sort_memory_per_thread,
                staging_dir=alignment_staging_dir,
            )
        )
        bam_path = aligner.align(entry)
        header_matches = aligner.verify_matches_reference_fasta(bam_path, args.reference_fasta)
    except Exception as e:
        notify_discord(f"MAT-20 install/alignment FAILED for sample {args.sample_title}: {e}")
        raise

    print(f"aligned BAM: {bam_path}")
    print(f"@SQ header matches {args.reference_fasta}: {header_matches}")
    if not header_matches:
        notify_discord(f"MAT-20 install/alignment FAILED for sample {args.sample_title}: @SQ header does not match {args.reference_fasta}")
        raise SystemExit(1)

    # Alignment already succeeded and was verified above -- a Discord hiccup
    # here must not discard that result (this is what silently killed IC17/IC37
    # after a full day of alignment: an unguarded call here raised and the
    # whole job died with no other trace).
    try:
        notify_discord(f"MAT-20 install/alignment finished for sample {args.sample_title}: {bam_path}")
    except Exception as e:
        print(f"WARNING: failed to send Discord notification: {e}")


if __name__ == "__main__":
    main()
