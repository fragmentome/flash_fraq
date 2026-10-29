from dataclasses import dataclass
from pathlib import Path

from pandas import DataFrame, Series, read_csv
from statsmodels.api import nonparametric

from flash_fraq.genome_mask import DEFAULT_MASK_MARGIN_BP, mask_frame

# DELFI window size (phase2 §3); used only as a fallback bin width for masking when a
# windows table lacks an explicit `end` column (e.g. hand-built test fixtures) --
# real `window_fragment_counts.tsv` rows always carry `end`.
DELFI_WINDOW_BP: int = 5_000_000

# Approximate hg19 centromere midpoints (Mb), used only to bucket 5 Mb windows into
# chromosome arms for the arm-level fragmentation z-score (phase2 §9 arm_zscore_*).
# Not a precise cytoband table -- fine for arm bucketing at 5 Mb window resolution,
# but should be validated against a real hg19 centromere/cytoband file before HPC use.
HG19_CENTROMERE_MB: dict[str, float] = {
    "chr1": 123.4, "chr2": 93.9, "chr3": 90.9, "chr4": 50.0, "chr5": 48.8,
    "chr6": 59.8, "chr7": 60.1, "chr8": 45.2, "chr9": 43.0, "chr10": 39.8,
    "chr11": 53.4, "chr12": 35.5, "chr13": 17.7, "chr14": 17.2, "chr15": 19.0,
    "chr16": 36.8, "chr17": 24.0, "chr18": 17.2, "chr19": 26.5, "chr20": 27.5,
    "chr21": 13.2, "chr22": 14.7, "chrX": 60.6,
}


def arm_id(chrom: str, start: int) -> str:
    """Chromosome-arm label ("chr1p"/"chr1q") for a 5 Mb window's start coordinate.

    Module-level so a caller can bucket windows into arms -- e.g. the arm_zscores.csv
    label back-fill in scripts/ -- without instantiating a corrector.
    """
    centromere_bp = HG19_CENTROMERE_MB[chrom] * 1_000_000
    return f"{chrom}{'p' if start < centromere_bp else 'q'}"


@dataclass(frozen=True)
class GcCorrectionConfig:
    input_windows_path: Path
    output_windows_path: Path
    window_ratios_path: Path
    arm_zscore_path: Path
    span: float = 0.75
    # Centromere/chromosome-end mask (docs/results-corruption-audit.md's coverage-bias
    # finding): windows within `mask_margin_bp` of either are dropped before the LOESS
    # fit and from every MAT-15 output. 0 disables masking entirely.
    mask_margin_bp: int = DEFAULT_MASK_MARGIN_BP

    @classmethod
    def from_sources(
        cls,
        input_windows_path: Path,
        output_windows_path: Path,
        window_ratios_path: Path | None = None,
        arm_zscore_path: Path | None = None,
        span: float = 0.75,
        mask_margin_bp: int = DEFAULT_MASK_MARGIN_BP,
    ) -> "GcCorrectionConfig":
        input_windows_path = Path(input_windows_path)
        output_windows_path = Path(output_windows_path)
        return cls(
            input_windows_path=input_windows_path,
            output_windows_path=output_windows_path,
            window_ratios_path=(
                Path(window_ratios_path)
                if window_ratios_path is not None
                else output_windows_path.parent / "window_ratios.csv"
            ),
            arm_zscore_path=(
                Path(arm_zscore_path)
                if arm_zscore_path is not None
                else output_windows_path.parent / "arm_zscores.csv"
            ),
            span=span,
            mask_margin_bp=mask_margin_bp,
        )


class GcCorrector:
    """Applies DELFI-style LOESS GC-bias correction to 5 Mb window fragment counts.

    Per phase2-preprocessing-plan.md §4: GC correction is fit and applied
    separately to short and long fragment counts (not directly to their
    ratio) via LOESS (span = 3/4) of count against window GC content — the
    predicted trend is subtracted from each count and the genome-wide
    median count is added back to keep counts on their original scale. The
    short:long ratio is then recomputed from the corrected counts.
    """

    def __init__(self, config: GcCorrectionConfig) -> None:
        self.config = config

    def run(self) -> DataFrame:
        windows = read_csv(self.config.input_windows_path, sep="\t")
        corrected = self.correct(windows)
        corrected.to_csv(self.config.output_windows_path, sep="\t", index=False)
        return corrected

    def run_full(self) -> tuple[Path, Path]:
        """Runs GC correction and writes both MAT-17-facing outputs: window_ratios.csv
        and arm_zscores.csv. Also writes the full diagnostic table to output_windows_path."""
        corrected = self.run()
        window_ratios_path = self.write_window_ratios(corrected)
        arm_zscore_path = self.write_arm_zscores(corrected)
        return window_ratios_path, arm_zscore_path

    def correct(self, windows: DataFrame) -> DataFrame:
        windows = self._drop_masked(windows)
        corrected = windows.copy()
        corrected["short_count_gc_corrected"] = self._correct_counts(
            windows["gc_content"], windows["short_count"]
        )
        corrected["long_count_gc_corrected"] = self._correct_counts(
            windows["gc_content"], windows["long_count"]
        )
        corrected["short_long_ratio"] = (
            corrected["short_count_gc_corrected"] / corrected["long_count_gc_corrected"]
        )
        corrected["gc_corrected_total_coverage"] = (
            corrected["short_count_gc_corrected"] + corrected["long_count_gc_corrected"]
        )
        corrected["gc_corrected_short_coverage"] = corrected["short_count_gc_corrected"]
        return corrected

    def write_window_ratios(self, corrected: DataFrame) -> Path:
        self.config.window_ratios_path.parent.mkdir(parents=True, exist_ok=True)
        window_ratios = corrected[
            ["sample_id", "window_index", "short_long_ratio", "gc_corrected_total_coverage", "gc_corrected_short_coverage"]
        ]
        window_ratios.to_csv(self.config.window_ratios_path, index=False)
        return self.config.window_ratios_path

    def compute_arm_zscores(self, corrected: DataFrame) -> DataFrame:
        """Chromosome-arm fragmentation z-scores (phase2 §9 arm_zscore_*).

        Each 5 Mb window is bucketed into a chromosome arm via its centromere
        position, the corrected short:long ratio is averaged per arm, and each
        arm's mean is z-scored against this sample's own across-arm distribution.

        Keyed on the `arm` label rather than on a positional index: a sample only has
        windows on the arms its coverage reaches (44, 45 or 46 of the 46 primary-contig
        arms across the real cohort), so a per-sample rank would put different arms in
        the same column of any cross-sample join.

        This is a self-normalizing simplification: DELFI's original arm z-scores
        are computed against a healthy-reference cohort's mean/std, not a sample's
        own arms. Revisit once a reference cohort baseline is available.
        """
        with_arm = corrected.copy()
        with_arm["arm"] = with_arm.apply(lambda row: self._arm_id(row["chrom"], row["start"]), axis=1)
        arm_means = with_arm.groupby("arm")["short_long_ratio"].mean().sort_index()
        arm_std = arm_means.std(ddof=0)
        if arm_std == 0:
            arm_zscores = arm_means * 0.0
        else:
            arm_zscores = (arm_means - arm_means.mean()) / arm_std

        sample_id = corrected["sample_id"].iloc[0]
        return DataFrame(
            {
                "sample_id": sample_id,
                "arm": arm_zscores.index,
                "zscore": arm_zscores.to_numpy(),
            }
        )

    def write_arm_zscores(self, corrected: DataFrame) -> Path:
        self.config.arm_zscore_path.parent.mkdir(parents=True, exist_ok=True)
        arm_zscores = self.compute_arm_zscores(corrected)
        arm_zscores[["sample_id", "arm", "zscore"]].to_csv(self.config.arm_zscore_path, index=False)
        return self.config.arm_zscore_path

    def _arm_id(self, chrom: str, start: int) -> str:
        return arm_id(chrom, start)

    def _drop_masked(self, windows: DataFrame) -> DataFrame:
        """Drops centromere-/chromosome-end-proximal windows (§ module docstring of
        `genome_mask`) before anything is fit or written.

        Because this runs first in `correct()`, the mask reaches every MAT-15 output
        for free: the LOESS fit (never sees a masked window's count), the diagnostic
        `output_windows_path` table, `window_ratios.csv` (`write_window_ratios` just
        selects columns off whatever `correct()` returned) and the arm z-score means
        (`compute_arm_zscores` groups over the same filtered frame). An arm on an
        acrocentric chromosome (13/14/15/21/22) whose whole p-arm falls inside the
        centromere margin loses every one of its windows here and so never appears in
        `arm_zscores.csv` -- the same "ragged arm set" a sample with low coverage on an
        arm already produces (see `compute_arm_zscores`'s docstring), not a new failure
        mode `feature_matrix.pivot_arm_zscores` needs to special-case.

        `end` may be absent from a hand-built windows table (e.g. a test fixture); the
        window is then assumed to be a full-width DELFI window starting at `start`.
        """
        if self.config.mask_margin_bp <= 0:
            return windows
        has_end = "end" in windows.columns
        is_masked = mask_frame(
            windows,
            chrom_col="chrom",
            start_col="start",
            end_col="end" if has_end else None,
            bin_size_bp=None if has_end else DELFI_WINDOW_BP,
            margin_bp=self.config.mask_margin_bp,
        )
        return windows.loc[~is_masked].reset_index(drop=True)

    def _correct_counts(self, gc_content: Series, counts: Series) -> Series:
        smoothed_trend = nonparametric.lowess(
            counts,
            gc_content,
            frac=self.config.span,
            return_sorted=False,
        )
        return counts - smoothed_trend + counts.median()
