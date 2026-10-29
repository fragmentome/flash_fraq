"""Generate the DELFI 5 Mb window and ichorCNA 1 Mb bin BEDs that
fragment_extraction.py takes as inputs (docs/phase2-preprocessing-plan.md §3).

Neither BED shipped with the repo -- only the test suite ever built them, from
synthetic contigs. Invoked by the localfeat_make_windows rule in
workflow/rules/local_features.smk.

Windows are named "<prefix><n>" with a genome-wide running index, because
FragmentExtractor._parse_window_index() reads the trailing integer off the
BED name column and uses it as window_index in window_fragment_counts.tsv.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path


def read_contig_lengths(fai_path: Path, contigs: list[str]) -> list[tuple[str, int]]:
    """Contig lengths in the order the caller listed them, not .fai order.

    hg19.fa.fai is ordered by size (chrX sits between chr7 and chr8), so
    trusting file order would interleave the windows and make window_index
    non-monotonic along the genome.
    """
    lengths: dict[str, int] = {}
    with fai_path.open() as handle:
        for line in handle:
            name, length = line.split("\t")[:2]
            lengths[name] = int(length)

    missing = [contig for contig in contigs if contig not in lengths]
    if missing:
        raise ValueError(f"contigs absent from {fai_path}: {', '.join(missing)}")
    return [(contig, lengths[contig]) for contig in contigs]


def write_windows(output_path: Path, contig_lengths: list[tuple[str, int]], window_bp: int, name_prefix: str) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    window_index = 0
    with output_path.open("w") as out:
        for contig, length in contig_lengths:
            for start in range(0, length, window_bp):
                # Trailing partial window is kept but clipped to the contig
                # end; GC content is computed over whatever bases it holds.
                end = min(start + window_bp, length)
                window_index += 1
                out.write(f"{contig}\t{start}\t{end}\t{name_prefix}{window_index}\n")
    return window_index


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--fai", required=True, type=Path)
    parser.add_argument("--contigs", required=True, nargs="+")
    parser.add_argument("--delfi-window-bp", required=True, type=int)
    parser.add_argument("--ichor-bin-bp", required=True, type=int)
    parser.add_argument("--delfi-bed", required=True, type=Path)
    parser.add_argument("--ichor-bed", required=True, type=Path)
    args = parser.parse_args()

    contig_lengths = read_contig_lengths(args.fai, args.contigs)

    delfi_count = write_windows(args.delfi_bed, contig_lengths, args.delfi_window_bp, "window")
    ichor_count = write_windows(args.ichor_bed, contig_lengths, args.ichor_bin_bp, "bin")

    print(f"DELFI windows ({args.delfi_window_bp} bp): {delfi_count} -> {args.delfi_bed}")
    print(f"ichor bins ({args.ichor_bin_bp} bp): {ichor_count} -> {args.ichor_bed}")


if __name__ == "__main__":
    main()
