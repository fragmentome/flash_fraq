"""Verify one staged CRAM before anything expensive reads it
(docs/phase2-preprocessing-plan.md §8b.3 step 1).

Two checks, for two different failure modes:

1. `samtools quickcheck` -- catches a truncated file. sync_crams_to_hdd.sh
   runs rsync with --partial and no --partial-dir, so an interrupted transfer
   leaves partial bytes at the *final* filename. Such a CRAM opens fine and
   decodes for a while before dying mid-stream, which without this gate would
   surface as a confusing failure hours into a fragment-extraction pass.

2. @SQ vs. the reference .fai -- closes the phase2 "reference FASTA variant"
   open item. A matching build *name* is not enough: HCC1395 has to be aligned
   against the identical contig set, or the `samtools merge` inside
   admixture.py's merge_and_dedup() will reject the mix.

Invoked by the localfeat_verify_cram rule in workflow/rules/local_features.smk.
"""

from __future__ import annotations

import json
import subprocess
from argparse import ArgumentParser
from pathlib import Path

import pysam


def quickcheck(cram_path: Path) -> None:
    result = subprocess.run(
        ["samtools", "quickcheck", "-v", str(cram_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools quickcheck rejected {cram_path} (exit {result.returncode}): "
            f"{result.stdout.strip()} {result.stderr.strip()} -- most likely an "
            f"incomplete rsync --partial transfer; re-run the sync for this sample."
        )


def read_cram_contigs(cram_path: Path, reference_fasta: Path) -> list[tuple[str, int]]:
    with pysam.AlignmentFile(str(cram_path), "rc", reference_filename=str(reference_fasta)) as cram:
        return list(zip(cram.references, cram.lengths))


def read_fai_contigs(fai_path: Path) -> list[tuple[str, int]]:
    contigs: list[tuple[str, int]] = []
    with fai_path.open() as handle:
        for line in handle:
            name, length = line.split("\t")[:2]
            contigs.append((name, int(length)))
    return contigs


def compare(cram_contigs: list[tuple[str, int]], fai_contigs: list[tuple[str, int]]) -> dict[str, object]:
    """Order matters, not just membership -- @SQ order defines the BAM/CRAM
    reference-id encoding, so two files with the same contigs in a different
    order still cannot be merged."""
    cram_names = [name for name, _ in cram_contigs]
    fai_names = [name for name, _ in fai_contigs]

    cram_lengths = dict(cram_contigs)
    fai_lengths = dict(fai_contigs)
    shared = set(cram_lengths) & set(fai_lengths)

    # cram_names == fai_names is only informative when the CRAM lists the whole
    # reference; a CRAM that legitimately carries a subset makes it False no
    # matter how the contigs are ordered, so on its own it is either trivially
    # true or uselessly false. order_matches_reference is the check that
    # actually means something: the CRAM's contigs, in the reference's order,
    # restricted to the contigs the CRAM has.
    reference_order_of_cram_contigs = [name for name in fai_names if name in set(cram_names)]

    return {
        "identical_order": cram_names == fai_names,
        "order_matches_reference": cram_names == reference_order_of_cram_contigs,
        "missing_from_reference": sorted(set(cram_names) - set(fai_names)),
        "missing_from_cram": sorted(set(fai_names) - set(cram_names)),
        "length_mismatches": [
            {"contig": name, "cram": cram_lengths[name], "reference": fai_lengths[name]}
            for name in sorted(shared)
            if cram_lengths[name] != fai_lengths[name]
        ],
    }


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--cram-path", required=True, type=Path)
    parser.add_argument("--reference-fai", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-path", required=True, type=Path)
    args = parser.parse_args()

    quickcheck(args.cram_path)

    cram_contigs = read_cram_contigs(args.cram_path, args.reference_fasta)
    fai_contigs = read_fai_contigs(args.reference_fai)
    comparison = compare(cram_contigs, fai_contigs)

    # A contig in the reference but not the CRAM is normal and harmless -- the
    # CRAM header only lists what the aligner emitted. The reverse is fatal:
    # reads placed on a contig the reference lacks cannot be decoded, and a
    # length mismatch means a different assembly patch under the same name.
    # So is a reordering: @SQ order defines the reference-id encoding, which is
    # the whole reason this check compares order rather than membership, and
    # recording a mismatch while still passing the sample would defeat it.
    fatal = (
        bool(comparison["missing_from_reference"])
        or bool(comparison["length_mismatches"])
        or not comparison["order_matches_reference"]
    )

    report = {
        "sample_id": args.sample_id,
        "cram_path": str(args.cram_path),
        "reference_fasta": str(args.reference_fasta),
        "quickcheck_passed": True,
        "cram_contig_count": len(cram_contigs),
        "reference_contig_count": len(fai_contigs),
        "passed": not fatal,
        **comparison,
    }

    args.output_path.parent.mkdir(parents=True, exist_ok=True)
    args.output_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))

    if fatal:
        raise SystemExit(
            f"{args.sample_id}: CRAM @SQ headers are incompatible with {args.reference_fasta}; "
            f"see {args.output_path}"
        )


if __name__ == "__main__":
    main()
