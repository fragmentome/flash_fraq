"""Assemble the Phase 3 output-contract feature matrix (docs/phase2-preprocessing-plan.md §9).

Joins the per-mixed-sample outputs of MAT-13 (empirical tumor fraction, sample/source/
platform metadata), MAT-14 (fragment/TLEN + window intersections), MAT-15 (GC-corrected
5 Mb window ratios), and MAT-16 (1 Mb CN bins + ichorCNA estimate) into a single
one-row-per-sample feature matrix.

`empirical_tumor_fraction` (not `nominal_tumor_fraction`) is the regression label
(phase2 §7.4). `ichorCNA_tumor_fraction_estimate` is ablation-only and excluded from
the primary feature set by default (phase0 §2, phase2 §9). `split_group` is assigned
at the source-subject level, before mixing (phase0 §5.2) -- by
`admixture_design.MixturePlanner`, which emits the `split_assignment.csv` read here.
"""

from dataclasses import dataclass
from pathlib import Path

from pandas import DataFrame, concat, merge, read_csv

from flash_fraq.genome_mask import DEFAULT_MASK_MARGIN_BP, mask_frame

# ichorCNA's configured bin size (cnv_binning.IchorCnaConfig.bin_size_bp's default);
# cnv_bins.csv carries each bin's `start` but not its `end`, so masking a bin needs this
# to reconstruct `end = start + CNV_BIN_BP` the same way genome_mask.mask_frame does for
# any other fixed-size-bin table.
CNV_BIN_BP: int = 1_000_000

ABLATION_ONLY_COLUMNS: tuple[str, ...] = ("ichorCNA_tumor_fraction_estimate",)

REGRESSION_LABEL_COLUMN: str = "empirical_tumor_fraction"

# Split group for a source subject that appears in *every* split by construction.
#
# phase0 §5.2 splits by source subject so no source leaks across splits. That is
# achievable for the healthy background pool (hundreds of subjects) but not for the
# tumor side: v1 has exactly one tumor source, HCC1395, and it is in every mixture.
# Assigning it to one split would confine the whole matrix to that split. So it is
# assigned this sentinel instead, and a mixed sample takes its split from whichever
# side is split-bound. The leakage guarantee §5.2 actually buys is therefore
# "no *healthy background* subject spans two splits", and holding out an entire
# source/platform (§5.2's second bullet) remains the check on the tumor side.
SPLIT_GROUP_ANY: str = "all"


def pivot_long_to_wide(
    long_df: DataFrame,
    column_key: str,
    value_key: str,
    column_prefix: str,
    index_key: str = "sample_id",
) -> DataFrame:
    """One-row-per-`index_key` pivot of a long table, with prefixed column names.

    Shared by `FeatureMatrixAssembler` (one row per mixed admixture sample) and by
    `cohort_matrix.CohortFeatureMatrixAssembler` (one row per real cohort sample), so
    both produce byte-identical column naming for the same upstream table. The index is
    left on the frame -- callers that want a mergeable table call `.reset_index()`.

    `pivot_table` rather than `pivot`: the cohort assembler concatenates hundreds of
    per-sample tables, and `pivot` raises on a duplicate (sample_id, key) pair instead
    of surfacing it as something the caller can check.
    """
    wide = long_df.pivot_table(index=index_key, columns=column_key, values=value_key, aggfunc="mean")
    wide.columns = [f"{column_prefix}{_format_key(key)}" for key in wide.columns]
    return wide


def _format_key(key: object) -> str:
    """Numeric pivot keys render as plain ints (`window_1`, not `window_1.0`); string
    keys (e.g. the cohort matrix's chromosome-arm labels) pass through unchanged."""
    if isinstance(key, str):
        return key
    return str(int(key))


def pivot_arm_zscores(long_df: DataFrame, column_key: str = "arm") -> DataFrame:
    """arm_zscores.csv -> `arm_zscore_*`, one column per arm (`arm_zscore_chr1p`).

    Keyed on the arm label, not on the arm's position in the file: the per-sample arm
    set is ragged, so a positional key would place different arms under one column.
    """
    return pivot_long_to_wide(long_df, column_key=column_key, value_key="zscore", column_prefix="arm_zscore_")


def pivot_window_ratios(long_df: DataFrame) -> DataFrame:
    """window_ratios.csv -> the three §9 per-window blocks, in §9's column order."""
    blocks = [
        pivot_long_to_wide(long_df, "window_index", "short_long_ratio", "short_long_ratio_window_"),
        pivot_long_to_wide(long_df, "window_index", "gc_corrected_total_coverage", "gc_corrected_total_coverage_window_"),
        pivot_long_to_wide(long_df, "window_index", "gc_corrected_short_coverage", "gc_corrected_short_coverage_window_"),
    ]
    return concat(blocks, axis=1, sort=False)


def pivot_cnv_bins(long_df: DataFrame, mask_margin_bp: int = DEFAULT_MASK_MARGIN_BP) -> DataFrame:
    """cnv_bins.csv -> `cnv_logr_bin_*`, dropping bins within `mask_margin_bp` of an hg19
    centromere or chromosome end (docs/results-corruption-audit.md's coverage-bias
    finding) before pivoting. 0 disables masking.

    Shared by `FeatureMatrixAssembler` (admixture) and `CohortFeatureMatrixAssembler`
    (real cohort) -- both read the same per-sample `cnv_bins.csv`/`ichor/cnv_bins.csv`
    written by `cnv_binning.IchorCnaRunner.write_cnv_bins`, which already carries `chrom`
    and `start` for exactly this reason (see its `_bin_indices` docstring), so both
    callers get the mask without a separate code path. The bin set masked out is the same
    for every sample (it depends only on genomic position), so it stays consistent across
    the columns any cross-sample join produces.
    """
    if mask_margin_bp > 0:
        is_masked = mask_frame(long_df, chrom_col="chrom", start_col="start", end_col=None, bin_size_bp=CNV_BIN_BP, margin_bp=mask_margin_bp)
        long_df = long_df.loc[~is_masked]
    return pivot_long_to_wide(long_df, column_key="bin_index", value_key="cnv_logr", column_prefix="cnv_logr_bin_")


@dataclass(frozen=True)
class FeatureMatrixConfig:
    sample_metadata_path: Path
    fragment_stats_path: Path
    arm_zscore_path: Path
    window_ratio_path: Path
    cnv_bin_path: Path
    ichor_estimate_path: Path
    split_assignment_path: Path
    output_path: Path
    # Centromere/chromosome-end mask applied to cnv_bins.csv before pivoting into
    # cnv_logr_bin_* (see pivot_cnv_bins). window_ratios.csv/arm_zscores.csv need no
    # separate knob here: gc_correction.GcCorrectionConfig.mask_margin_bp already dropped
    # those windows upstream, before those files were written.
    mask_margin_bp: int = DEFAULT_MASK_MARGIN_BP

    @classmethod
    def from_sources(
        cls, input_dir: Path, output_path: Path, mask_margin_bp: int = DEFAULT_MASK_MARGIN_BP
    ) -> "FeatureMatrixConfig":
        return cls(
            sample_metadata_path=input_dir / "sample_metadata.csv",
            fragment_stats_path=input_dir / "fragment_stats.csv",
            arm_zscore_path=input_dir / "arm_zscores.csv",
            window_ratio_path=input_dir / "window_ratios.csv",
            cnv_bin_path=input_dir / "cnv_bins.csv",
            ichor_estimate_path=input_dir / "ichor_estimates.csv",
            split_assignment_path=input_dir / "split_assignment.csv",
            output_path=output_path,
            mask_margin_bp=mask_margin_bp,
        )


class FeatureMatrixAssembler:
    """Builds the one-row-per-mixed-sample §9 feature matrix from upstream tables."""

    def __init__(self, config: FeatureMatrixConfig) -> None:
        self.config = config

    def load_sample_metadata(self) -> DataFrame:
        return read_csv(self.config.sample_metadata_path)

    def load_fragment_stats(self) -> DataFrame:
        return read_csv(self.config.fragment_stats_path)

    def load_arm_zscores_wide(self) -> DataFrame:
        return pivot_arm_zscores(read_csv(self.config.arm_zscore_path)).reset_index()

    def load_window_ratios_wide(self) -> DataFrame:
        return pivot_window_ratios(read_csv(self.config.window_ratio_path)).reset_index()

    def load_cnv_bins_wide(self) -> DataFrame:
        return pivot_cnv_bins(read_csv(self.config.cnv_bin_path), mask_margin_bp=self.config.mask_margin_bp).reset_index()

    def load_ichor_estimates(self) -> DataFrame:
        return read_csv(self.config.ichor_estimate_path)

    def load_split_assignment(self) -> DataFrame:
        return read_csv(self.config.split_assignment_path)

    def assign_split_groups(self, sample_metadata: DataFrame, split_assignment: DataFrame) -> DataFrame:
        """Resolves per-sample `split_group` from a subject-level assignment made before mixing.

        A mixed sample's tumor and healthy source subjects must already share a split
        group (phase0 §5.2) — a mismatch means a source subject leaked across splits.
        """
        subject_to_split = split_assignment.set_index("subject_id")["split_group"]

        tumor_split = sample_metadata["tumor_subject_id"].map(subject_to_split)
        healthy_split = sample_metadata["healthy_subject_id"].map(subject_to_split)

        # A side assigned SPLIT_GROUP_ANY defers to the other one; where both sides are
        # split-bound they must still agree, which is the original leakage check.
        resolved_tumor = tumor_split.where(tumor_split != SPLIT_GROUP_ANY, healthy_split)
        resolved_healthy = healthy_split.where(healthy_split != SPLIT_GROUP_ANY, tumor_split)

        if not resolved_tumor.equals(resolved_healthy):
            mismatched = sample_metadata.loc[resolved_tumor != resolved_healthy, "sample_id"]
            raise ValueError(
                f"tumor/healthy source subjects assigned to different split groups for samples: "
                f"{mismatched.tolist()}"
            )

        result = sample_metadata.copy()
        result["split_group"] = resolved_tumor
        return result

    def assemble(self) -> DataFrame:
        sample_metadata = self.load_sample_metadata()
        split_assignment = self.load_split_assignment()
        sample_metadata = self.assign_split_groups(sample_metadata, split_assignment)

        matrix = sample_metadata
        for other in (
            self.load_fragment_stats(),
            self.load_arm_zscores_wide(),
            self.load_window_ratios_wide(),
            self.load_cnv_bins_wide(),
            self.load_ichor_estimates(),
        ):
            matrix = merge(matrix, other, on="sample_id", how="left")

        return matrix

    def write(self, matrix: DataFrame) -> Path:
        self.config.output_path.parent.mkdir(parents=True, exist_ok=True)
        matrix.to_csv(self.config.output_path, index=False)
        return self.config.output_path

    def get_primary_feature_columns(self, matrix: DataFrame) -> list[str]:
        """Feature columns for the primary model — excludes ablation-only columns."""
        non_feature_columns = {
            "sample_id",
            "tumor_source_id",
            "healthy_source_id",
            "platform_tumor",
            "platform_healthy",
            "reference_build",
            "target_coverage",
            "tumor_subject_id",
            "healthy_subject_id",
            "nominal_tumor_fraction",
            REGRESSION_LABEL_COLUMN,
            # See tf_model.NON_FEATURE_COLUMNS: empirical_tumor_fraction_sites *is*
            # REGRESSION_LABEL_COLUMN's value verbatim in this cohort (the canonical
            # cross-estimator label), and empirical_tumor_fraction_clonal is the other
            # estimator on the same SNVs -- both leak the label, never features.
            "empirical_tumor_fraction_sites",
            "empirical_tumor_fraction_clonal",
            "split_group",
            # audit finding 2: HCC1395 gDNA carries ~10^-2.6 mtDNA fraction vs cfDNA's
            # ~10^-4.9, so this is a tumor-read counter (rho=0.96 vs nominal f in
            # real_background), not a fragmentomics measurement. Stays in the matrix
            # as metadata, never a candidate feature.
            "mtdna_fraction_log10",
        }
        excluded = non_feature_columns | set(ABLATION_ONLY_COLUMNS)
        return [column for column in matrix.columns if column not in excluded]
