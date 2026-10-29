"""Download one HCC1395/HCC1395BL run's FASTQ pair -- from ENA's submitted
fastq.gz where it's served that way, else the SRA toolkit
(docs/phase1-sourcing.md; MAT-23).

Invoked by workflow/rules/hcc1395.smk's `fetch_hcc1395_fastq` rule.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.hcc1395_retrieval import HCC1395FastqRetriever, HCC1395RetrievalConfig
from flash_fraq.notifications import notify_discord


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--run-accession", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--threads", type=int, default=4)
    # Escape hatch for forcing the SRA toolkit path on a run ENA does serve
    # (e.g. to compare the two); availability is otherwise checked per run.
    parser.add_argument("--no-ena", dest="prefer_ena", action="store_false")
    args = parser.parse_args()

    retriever = HCC1395FastqRetriever(
        HCC1395RetrievalConfig.from_sources(
            output_dir=args.output_dir, threads=args.threads, prefer_ena=args.prefer_ena
        )
    )
    try:
        record = retriever.download_run(args.run_accession)
    except Exception as e:
        notify_discord(f"MAT-23 HCC1395 fetch FAILED for run {args.run_accession}: {e}")
        raise

    print(f"retrieval source: {record.source}")
    print(f"fastq R1: {record.r1_path}")
    print(f"fastq R2: {record.r2_path}")


if __name__ == "__main__":
    main()
