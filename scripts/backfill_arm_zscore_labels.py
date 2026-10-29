"""One-off migration: add the `arm` label column to already-written arm_zscores.csv.

`GcCorrector.write_arm_zscores` used to drop the `arm` column it had already
computed and write only (sample_id, arm_index, zscore). `arm_index` is a dense
per-sample rank over the arms *that sample* has windows on, so it is not a
cross-sample identifier -- see docs/todo/in-progress/arm-label-writer-fix-and-backfill.md.
The writer is fixed; this back-fills the samples produced before the fix.

The label is recomputed from each sample's own window_fragment_counts.gc_corrected.tsv
via `gc_correction.arm_id`, which is exactly what the writer bucketed on, so no GC
correction is re-run. A sample is rewritten only if the recomputed arm set has the
same size as its existing arm_zscores.csv; a disagreement means the two files came
from different window data, and rewriting would invent labels rather than recover them.

Idempotent: a file that already carries `arm` is left alone, so the script can be
re-run over a partially migrated tree. Every existing field is carried through as
*text* (`dtype=str`), so the migration adds a column and changes nothing else: pandas'
CSV float parser is not correctly rounded, so parsing and re-writing a z-score would
shift a handful of them by an ulp and make the cohort matrix diff non-empty for reasons
that have nothing to do with arm labels.

  python scripts/backfill_arm_zscore_labels.py --derived-root /mnt/8tb_hdd/flash_fraq_derived \
      --dataset Snyder_Shendure_WPS --dataset Cristiano_Velculescu_DELFI
"""

from __future__ import annotations

import sys
from argparse import ArgumentParser
from dataclasses import dataclass
from pathlib import Path

from pandas import DataFrame, read_csv

from flash_fraq.gc_correction import arm_id

BACKFILLED_COLUMNS: tuple[str, ...] = ("sample_id", "arm_index", "arm", "zscore")


@dataclass(frozen=True)
class BackfillOutcome:
    sample_id: str
    status: str  # "backfilled" | "already_labelled" | "skipped" | "mismatch"
    detail: str = ""


class ArmLabelBackfiller:
    """Rewrites one derived tree's arm_zscores.csv files with their `arm` label."""

    def __init__(self, derived_root: Path, datasets: tuple[str, ...], dry_run: bool = False) -> None:
        self.derived_root = derived_root
        self.datasets = datasets
        self.dry_run = dry_run

    def sample_dirs(self) -> list[Path]:
        directories: list[Path] = []
        for dataset in self.datasets:
            dataset_dir = self.derived_root / dataset
            if not dataset_dir.is_dir():
                raise ValueError(f"dataset directory does not exist: {dataset_dir}")
            directories.extend(sorted(path for path in dataset_dir.iterdir() if path.is_dir()))
        return directories

    def arms_present(self, gc_corrected_windows_path: Path) -> list[str]:
        """The sample's arm set in the order `compute_arm_zscores` numbered it.

        That order is `groupby("arm").mean().sort_index()`, i.e. lexicographic over the
        arms present -- so `sorted()` here reproduces arm_index 1..N exactly.
        """
        windows = read_csv(gc_corrected_windows_path, sep="\t", usecols=["chrom", "start"])
        return sorted({arm_id(chrom, start) for chrom, start in zip(windows["chrom"], windows["start"])})

    def backfill_sample(self, sample_dir: Path) -> BackfillOutcome:
        sample_id = sample_dir.name
        arm_zscore_path = sample_dir / "arm_zscores.csv"
        windows_path = sample_dir / "window_fragment_counts.gc_corrected.tsv"
        if not arm_zscore_path.exists():
            return BackfillOutcome(sample_id, "skipped", "no arm_zscores.csv")

        zscores = read_csv(arm_zscore_path, dtype=str)
        if "arm" in zscores.columns:
            return BackfillOutcome(sample_id, "already_labelled")
        if not windows_path.exists():
            return BackfillOutcome(sample_id, "skipped", "no window_fragment_counts.gc_corrected.tsv")

        arms = self.arms_present(windows_path)
        if len(arms) != len(zscores):
            return BackfillOutcome(
                sample_id,
                "mismatch",
                f"arm_zscores.csv has {len(zscores)} arms but its GC-corrected window table implies {len(arms)}",
            )
        expected_index = [str(index) for index in range(1, len(arms) + 1)]
        if list(zscores["arm_index"]) != expected_index:
            return BackfillOutcome(
                sample_id,
                "mismatch",
                f"arm_index is not the dense rank 1..{len(arms)} the writer emits",
            )

        labelled = zscores.assign(arm=arms)[list(BACKFILLED_COLUMNS)]
        if not self.dry_run:
            self._write_atomically(labelled, arm_zscore_path)
        return BackfillOutcome(sample_id, "backfilled", f"{len(arms)} arms")

    def run(self) -> list[BackfillOutcome]:
        return [self.backfill_sample(sample_dir) for sample_dir in self.sample_dirs()]

    @staticmethod
    def _write_atomically(labelled: DataFrame, arm_zscore_path: Path) -> None:
        """Write via a sibling temp file + rename so an interrupted run never leaves a
        half-written table behind -- this rewrites the pipeline's own output in place."""
        temp_path = arm_zscore_path.with_suffix(".csv.backfill-tmp")
        labelled.to_csv(temp_path, index=False)
        temp_path.replace(arm_zscore_path)


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--derived-root", required=True, type=Path)
    parser.add_argument("--dataset", required=True, action="append", metavar="NAME")
    parser.add_argument("--dry-run", action="store_true", help="report what would change, write nothing")
    args = parser.parse_args()

    backfiller = ArmLabelBackfiller(args.derived_root, tuple(args.dataset), dry_run=args.dry_run)
    outcomes = backfiller.run()

    counts: dict[str, int] = {}
    for outcome in outcomes:
        counts[outcome.status] = counts.get(outcome.status, 0) + 1
    for status in ("backfilled", "already_labelled", "skipped", "mismatch"):
        print(f"{status}: {counts.get(status, 0)}")

    problems = [outcome for outcome in outcomes if outcome.status in {"skipped", "mismatch"}]
    for outcome in problems:
        print(f"  {outcome.status} {outcome.sample_id}: {outcome.detail}", file=sys.stderr)
    if any(outcome.status == "mismatch" for outcome in outcomes):
        raise SystemExit("refusing to complete: at least one sample's arm count disagrees with its window table")


if __name__ == "__main__":
    main()
