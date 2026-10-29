"""Driver for MAT-15 LOESS GC correction over one sample's 5 Mb window counts
(docs/phase2-preprocessing-plan.md §8b.3 step 6, DELFI half).

Operates purely on the small per-window table MAT-14 wrote -- it never opens
the alignment, so it is the one local stage with no staging or IO concern.

Invoked per-sample by the localfeat_gc_correct rule in
workflow/rules/local_features.smk.
"""

from __future__ import annotations

from argparse import ArgumentParser
from pathlib import Path

from flash_fraq.fragment_extraction import cleanup_fragment_intermediates
from flash_fraq.gc_correction import GcCorrectionConfig, GcCorrector
from flash_fraq.genome_mask import DEFAULT_MASK_MARGIN_BP


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--input-windows-path", required=True, type=Path)
    parser.add_argument("--output-windows-path", required=True, type=Path)
    parser.add_argument("--window-ratios-path", required=True, type=Path)
    parser.add_argument("--arm-zscore-path", required=True, type=Path)
    parser.add_argument("--span", type=float, default=0.75)
    parser.add_argument(
        "--mask-margin-bp",
        type=int,
        default=DEFAULT_MASK_MARGIN_BP,
        help="Drop windows within this many bp of a centromere or chromosome end before "
        "the LOESS fit (0 disables masking).",
    )
    parser.add_argument(
        "--cleanup-fragment-beds",
        action="store_true",
        help="After a successful run, remove fragment_midpoints.bed/fragments_delfi_5mb.bed/"
        "fragments_ichor_1mb.bed from --input-windows-path's directory (MAT-14's per-fragment "
        "intermediates -- see flash_fraq.fragment_extraction.cleanup_fragment_intermediates). "
        "Off by default; the localfeat_gc_correct Snakemake rule never passes it because its "
        "upstream localfeat_fragments rule already marks those files temp() and Snakemake "
        "deletes them itself. Opt in for drivers (e.g. admixture_mixture_gc_correct_driver.sh) "
        "that call this script directly instead of through Snakemake.",
    )
    args = parser.parse_args()

    corrector = GcCorrector(
        GcCorrectionConfig.from_sources(
            input_windows_path=args.input_windows_path,
            output_windows_path=args.output_windows_path,
            window_ratios_path=args.window_ratios_path,
            arm_zscore_path=args.arm_zscore_path,
            span=args.span,
            mask_margin_bp=args.mask_margin_bp,
        )
    )
    window_ratios_path, arm_zscore_path = corrector.run_full()

    print(f"window ratios: {window_ratios_path}")
    print(f"arm zscores: {arm_zscore_path}")

    if args.cleanup_fragment_beds:
        # window_ratios_path/arm_zscore_path just came back from a successful run_full()
        # above, so GC correction has already read and consumed the window-counts table
        # that was the .bed files' last purpose. Safe per
        # flash_fraq.fragment_extraction.FRAGMENT_BED_INTERMEDIATES.
        removed = cleanup_fragment_intermediates(args.input_windows_path.parent)
        for path in removed:
            print(f"removed fragment-bed intermediate: {path}")


if __name__ == "__main__":
    main()
