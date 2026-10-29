"""Driver for MAT-14 fragment extraction over one staged CRAM
(docs/phase2-preprocessing-plan.md §8b.3 step 5).

Runs on the workstation against the SSD-staged copy, not on the cluster and
not against the HDD mirror -- see the module docstring in
workflow/rules/local_features.smk for why staging comes first.

Invoked per-sample by the localfeat_fragments rule.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.fragment_extraction import FragmentExtractionConfig, FragmentExtractor


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--input-alignment", required=True, type=Path, help="staged .cram (or .bam)")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--delfi-windows-bed", required=True, type=Path)
    parser.add_argument("--ichor-bins-bed", required=True, type=Path)
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    args = parser.parse_args()

    extractor = FragmentExtractor(
        FragmentExtractionConfig.from_sources(
            input_bam=args.input_alignment,
            output_dir=args.output_dir,
            delfi_windows_bed=args.delfi_windows_bed,
            ichor_bins_bed=args.ichor_bins_bed,
            reference_fasta=args.reference_fasta,
            sample_id=args.sample_id,
        )
    )
    delfi_bed, ichor_bed, window_counts, fragment_stats = extractor.run_full()

    print(f"delfi intersected: {delfi_bed}")
    print(f"ichor intersected: {ichor_bed}")
    print(f"window fragment counts: {window_counts}")
    print(f"fragment stats: {fragment_stats}")


if __name__ == "__main__":
    main()
