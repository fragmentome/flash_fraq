"""Post-hoc tumor-read filter on existing admixture mixtures
(docs/todo/todo/reshaping-pilot-post-hoc-filter.md).

`fragment_reshaping.reshape_bam`'s accept/reject decision depends only on
`(seed, qname, tlen)`, an independent per-pair thinning -- so thinning the
tumor reads of an *already-built* mixture's `dedup.bam` after the fact is
equivalent in distribution to reshaping the tumor source before mixing, and
skips rebuilding the ~2.2-day mixture set. Each mixture's `dedup.bam` keeps
the tumor's and the background's own distinct `@RG` `SM` tags
(`admixture.py`'s `compute_empirical_tumor_fraction` docstring), so
`fragment_reshaping.tumor_rg_ids` tells the two apart within one merged BAM.

Two known side effects, both handled here per the ticket:
- Tumor depth alone would drop by the source-wide mean acceptance (more at
  high `f`), leaking `f` through the total-depth features. Fix: thin the
  background reads of every mixture in the batch by a per-mixture `k` so all
  of them land on one common total-pair target -- the most-thinned ("worst
  case") mixture's own post-thinning total.
- The realised admixture fraction becomes `f*a / (f*a + (1-f))` where `a` is
  this mixture's own realised acceptance -- the empirical SNV-VAF label
  already captures this; nothing to do here.

Runs in two passes per batch, both driven by `fragment_reshaping`:
1. `count_mixture_reshape_pairs` on every mixture (cheap, no output BAM) to
   size the common background-thinning target across the whole batch.
2. `reshape_mixture_bam` on every mixture, writing
   `<output_dir>/<sample_name>.reshaped.dedup.bam` + its index, and a summary
   CSV recording tumor/background counts, `ā_realised`, and `k` per mixture.

Never touches or overwrites an original `<sample_name>.dedup.bam` -- the
output name and (by construction, see `main`) directory always differ from
the input's.
"""

from __future__ import annotations

import csv
import json
from argparse import ArgumentParser
from pathlib import Path

import pysam

from flash_fraq.fragment_reshaping import (
    DEFAULT_MAX_INSERT_SIZE,
    DEFAULT_MIN_INSERT_SIZE,
    DEFAULT_MIN_MAPQ,
    compute_background_keep_threshold,
    count_mixture_reshape_pairs,
    read_acceptance_table_csv,
    reshape_mixture_bam,
    tumor_rg_ids as compute_tumor_rg_ids,
)

SUMMARY_CSV_COLUMNS: tuple[str, ...] = (
    "sample_name",
    "input_bam",
    "output_bam",
    "tumor_original",
    "tumor_kept",
    "background_original",
    "background_kept",
    "realised_acceptance",
    "background_keep_k",
    "total_kept",
    "common_target_total",
)


def sample_name_from_dedup_path(path: Path) -> str:
    """`mix_X_f0p500.dedup.bam` -> `mix_X_f0p500` -- strips the conventional
    `admixture.py` dedup suffix; falls back to the bare stem for any other
    naming so the script still works against a non-standard input path."""
    if path.name.endswith(".dedup.bam"):
        return path.name[: -len(".dedup.bam")]
    return path.stem


def load_mixture_paths(args) -> list[Path]:
    paths = [Path(p) for p in (args.mixture_bam or [])]
    if args.mixture_list is not None:
        with Path(args.mixture_list).open() as handle:
            paths.extend(Path(line.strip()) for line in handle if line.strip() and not line.startswith("#"))
    if not paths:
        raise ValueError("no mixture BAMs given: pass --mixture-bam (repeatable) and/or --mixture-list")
    return paths


def resolve_tumor_rg_ids(bam_path: Path, reference_fasta: Path, tumor_sample_name: str) -> set[str]:
    """An empty result is legitimate, not an error: an f=0 mixture is 100%
    background/normal reads (e.g. a matched_pair rung whose @RG SM is the
    normal sample, not the tumor), so `tumor_rg_ids` returning `set()` for it
    is expected -- `count_mixture_reshape_pairs`/`reshape_mixture_bam` treat
    every read not in `tumor_rg_ids` as background, so nothing is thinned at
    the tumor stage and the mixture just goes through the common background
    normalisation like every other mixture in the batch."""
    with pysam.AlignmentFile(str(bam_path), "rb", reference_filename=str(reference_fasta)) as bam:
        rg_ids = compute_tumor_rg_ids(bam.header, tumor_sample_name)
    if not rg_ids:
        print(
            f"warning: {bam_path}: no @RG with SM=={tumor_sample_name!r} found in the header "
            "-- treating all reads as background (expected for an f=0 mixture)"
        )
    return rg_ids


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mixture-bam", action="append", type=Path, help="One mixture's dedup.bam. Repeatable."
    )
    parser.add_argument(
        "--mixture-list", type=Path, help="Text file, one mixture dedup.bam path per line, combined with --mixture-bam."
    )
    parser.add_argument("--tumor-sample-name", required=True, help="The tumor source's @RG SM tag, e.g. HCC1395_WGS_NV_3.")
    parser.add_argument("--acceptance-table", required=True, type=Path, help="CSV from write_acceptance_table_csv (insert_size,probability).")
    parser.add_argument("--reference-fasta", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path, help="Must differ from every input BAM's directory.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-insert-size", type=int, default=DEFAULT_MIN_INSERT_SIZE)
    parser.add_argument("--max-insert-size", type=int, default=DEFAULT_MAX_INSERT_SIZE)
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=DEFAULT_MIN_MAPQ,
        help="Minimum MAPQ (both mates) for a tumor-RG pair to be retained; background-RG pairs are "
        "never MAPQ-filtered. See coverage_bias_report.json.",
    )
    parser.add_argument("--summary-csv", type=Path, default=None, help="Default: <output-dir>/reshape_admixture_mixture_summary.csv")
    args = parser.parse_args()

    mixture_paths = load_mixture_paths(args)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    for bam_path in mixture_paths:
        if bam_path.resolve().parent == output_dir.resolve():
            raise ValueError(
                f"--output-dir {output_dir} matches {bam_path}'s own directory; refusing to run "
                "anywhere near an original dedup.bam"
            )

    acceptance_table = read_acceptance_table_csv(args.acceptance_table)

    # Pass 1: count only, to size the common background-thinning target.
    per_mixture_rg_ids: dict[Path, set[str]] = {}
    counts_by_path = {}
    for bam_path in mixture_paths:
        rg_ids = resolve_tumor_rg_ids(bam_path, args.reference_fasta, args.tumor_sample_name)
        per_mixture_rg_ids[bam_path] = rg_ids
        counts_by_path[bam_path] = count_mixture_reshape_pairs(
            bam_path,
            args.reference_fasta,
            rg_ids,
            acceptance_table,
            args.seed,
            args.min_insert_size,
            args.max_insert_size,
            args.min_mapq,
        )

    totals_after_tumor_thin = {
        bam_path: counts.tumor_kept + counts.background_original for bam_path, counts in counts_by_path.items()
    }
    common_target_total = min(totals_after_tumor_thin.values())

    # Pass 2: write the reshaped BAM for each mixture, thinning background by
    # the per-mixture k that brings it down to common_target_total.
    rows = []
    for bam_path in mixture_paths:
        counts = counts_by_path[bam_path]
        needed_background = max(0, common_target_total - counts.tumor_kept)
        k = needed_background / counts.background_original if counts.background_original > 0 else 0.0
        threshold = compute_background_keep_threshold(
            bam_path, args.reference_fasta, per_mixture_rg_ids[bam_path], needed_background, args.seed
        )

        sample_name = sample_name_from_dedup_path(bam_path)
        output_bam = output_dir / f"{sample_name}.reshaped.dedup.bam"
        _, final_counts = reshape_mixture_bam(
            bam_path,
            args.reference_fasta,
            output_bam,
            per_mixture_rg_ids[bam_path],
            acceptance_table,
            threshold,
            args.seed,
            args.min_insert_size,
            args.max_insert_size,
            args.min_mapq,
        )
        pysam.index(str(output_bam))

        rows.append(
            {
                "sample_name": sample_name,
                "input_bam": str(bam_path),
                "output_bam": str(output_bam),
                "tumor_original": final_counts.tumor_original,
                "tumor_kept": final_counts.tumor_kept,
                "background_original": final_counts.background_original,
                "background_kept": final_counts.background_kept,
                "realised_acceptance": final_counts.realised_acceptance,
                "background_keep_k": k,
                "total_kept": final_counts.total_kept,
                "common_target_total": common_target_total,
            }
        )

    summary_csv = args.summary_csv or (output_dir / "reshape_admixture_mixture_summary.csv")
    with summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(json.dumps({"summary_csv": str(summary_csv), "n_mixtures": len(rows), "common_target_total": common_target_total, "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
