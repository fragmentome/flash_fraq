"""Cohort feature matrix for the *real* Snyder + Cristiano samples (pre-admixture).

`feature_matrix.py` assembles the phase2 §9 output contract: one row per **mixed
admixture sample**, keyed on `empirical_tumor_fraction` and carrying a
tumor/healthy source pair. That matrix cannot exist yet -- admixture needs
HCC1395 as the tumor source, which is not aligned.

What *does* exist is the per-sample MAT-14/15/16 output that
`workflow/rules/local_features.smk` already wrote for every mirrored CRAM, under
`<derived_root>/<dataset>/<sample_id>/`. This module pivots and joins that into one
cohort table, one row per **real single-source sample**, tagged with the dataset it
came from. Deliberate differences from the §9 matrix:

- No `empirical_tumor_fraction`: there is no admixture and no SNV-VAF ground truth.
  `ichorCNA_tumor_fraction_estimate` stands in wherever a fraction is needed (the §6
  correlation check), which makes that check a *weaker* proxy -- see
  `decoy_qc.DecoyClassifierQcConfig.for_cohort`.
- One `source_dataset`/`platform` pair per row instead of the tumor/healthy pair, since
  every row comes from a single library on a single platform.
- Arm z-scores are keyed by **arm label**, not by position: the per-sample arm set is
  ragged, so only the `arm` column in arm_zscores.csv is comparable across samples.

The pivot/column-naming helpers come from `feature_matrix`, so a column that survives
into Phase 3's real §9 matrix is named identically here.
"""

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

from pandas import DataFrame, concat, read_csv

from flash_fraq.feature_matrix import pivot_arm_zscores, pivot_cnv_bins, pivot_window_ratios

COHORT_TUMOR_FRACTION_COLUMN: str = "ichorCNA_tumor_fraction_estimate"

COHORT_METADATA_COLUMNS: tuple[str, ...] = ("sample_id", "source_dataset", "platform", "reference_build")

FRAGMENT_STATS_COLUMNS: tuple[str, ...] = (
    "fragment_size_median",
    "fragment_size_mode",
    "fragment_size_entropy",
)


@dataclass(frozen=True)
class CohortDataset:
    """One derived-output dataset directory and the platform label its samples carry.

    `platform` is metadata about the source, not something recoverable from the derived
    tables -- it is the §6 decoy classifier's target, so it is configured explicitly
    rather than inferred.
    """

    name: str
    platform: str


@dataclass(frozen=True)
class CohortSample:
    dataset: str
    sample_id: str
    platform: str
    sample_dir: Path

    @property
    def fragment_stats_path(self) -> Path:
        return self.sample_dir / "fragment_stats.csv"

    @property
    def window_ratio_path(self) -> Path:
        return self.sample_dir / "window_ratios.csv"

    @property
    def arm_zscore_path(self) -> Path:
        return self.sample_dir / "arm_zscores.csv"

    @property
    def gc_corrected_windows_path(self) -> Path:
        return self.sample_dir / "window_fragment_counts.gc_corrected.tsv"

    @property
    def cnv_bin_path(self) -> Path:
        return self.sample_dir / "ichor" / "cnv_bins.csv"

    @property
    def ichor_estimate_path(self) -> Path:
        return self.sample_dir / "ichor" / "ichor_estimates.csv"

    @property
    def complete_flag_path(self) -> Path:
        return self.sample_dir / "localfeat_sample_complete.flag"

    def required_paths(self) -> tuple[Path, ...]:
        return (
            self.fragment_stats_path,
            self.window_ratio_path,
            self.arm_zscore_path,
            self.gc_corrected_windows_path,
            self.cnv_bin_path,
            self.ichor_estimate_path,
        )


@dataclass(frozen=True)
class CohortMatrixConfig:
    derived_root: Path
    datasets: tuple[CohortDataset, ...]
    output_path: Path
    reference_build: str = "hg19"
    require_complete_flag: bool = True

    @classmethod
    def from_sources(
        cls,
        derived_root: Path,
        dataset_platforms: dict[str, str],
        output_path: Path,
        reference_build: str = "hg19",
    ) -> "CohortMatrixConfig":
        return cls(
            derived_root=Path(derived_root),
            datasets=tuple(CohortDataset(name=name, platform=platform) for name, platform in dataset_platforms.items()),
            output_path=Path(output_path),
            reference_build=reference_build,
        )


class CohortFeatureMatrixAssembler:
    """Builds the one-row-per-real-sample cohort matrix from the local_features outputs."""

    def __init__(self, config: CohortMatrixConfig) -> None:
        self.config = config

    def discover_samples(self) -> list[CohortSample]:
        """Every sample directory whose local_features run finished.

        Gated on `localfeat_sample_complete.flag` for the same reason
        `common.smk:discover_local_crams` gates on the .crai: the derived tree is
        written incrementally, and a sample whose ichor job is still running has a
        readable but partial set of tables. The per-file existence check below then
        catches a flag left behind by an older, narrower pipeline version.
        """
        discovered: list[CohortSample] = []
        seen: dict[str, str] = {}
        for dataset in self.config.datasets:
            dataset_dir = self.config.derived_root / dataset.name
            if not dataset_dir.is_dir():
                continue
            for sample_dir in sorted(path for path in dataset_dir.iterdir() if path.is_dir()):
                sample = CohortSample(
                    dataset=dataset.name,
                    sample_id=sample_dir.name,
                    platform=dataset.platform,
                    sample_dir=sample_dir,
                )
                if self.config.require_complete_flag and not sample.complete_flag_path.exists():
                    continue
                if not all(path.exists() for path in sample.required_paths()):
                    continue
                if sample.sample_id in seen:
                    raise ValueError(
                        f"sample_id {sample.sample_id!r} appears in both {seen[sample.sample_id]!r} and "
                        f"{dataset.name!r}; sample_id is the cohort matrix's join key and must be unique"
                    )
                seen[sample.sample_id] = dataset.name
                discovered.append(sample)
        return discovered

    def load_sample_metadata(self, samples: Iterable[CohortSample]) -> DataFrame:
        return DataFrame(
            [
                {
                    "sample_id": sample.sample_id,
                    "source_dataset": sample.dataset,
                    "platform": sample.platform,
                    "reference_build": self.config.reference_build,
                }
                for sample in samples
            ],
            columns=list(COHORT_METADATA_COLUMNS),
        )

    def load_long_table(self, samples: Iterable[CohortSample], path_attribute: str) -> DataFrame:
        """Concatenate one per-sample CSV across the cohort into a single long table."""
        frames = [read_csv(getattr(sample, path_attribute)) for sample in samples]
        if not frames:
            return DataFrame(columns=["sample_id"])
        return concat(frames, ignore_index=True)

    def load_arm_zscores_long(self, samples: Iterable[CohortSample]) -> DataFrame:
        """Per-sample arm z-scores, keyed on the stable `arm` label the writer records.

        `arm` is what makes a cohort-wide pivot meaningful: samples carry 44, 45 or 46 of
        the 46 primary-contig arms, so a row's position in the file is a dense per-sample
        rank rather than an identifier. Files written before the label existed carry that
        rank as `arm_index` and are rejected here rather than pivoted on -- the fix is
        scripts/backfill_arm_zscore_labels.py, not a fallback, so that a label can only
        ever come from one place.
        """
        frames: list[DataFrame] = []
        for sample in samples:
            zscores = read_csv(sample.arm_zscore_path)
            if "arm" not in zscores.columns:
                raise ValueError(
                    f"{sample.sample_id}: arm_zscores.csv predates the `arm` label column and cannot be "
                    f"keyed across samples; run scripts/backfill_arm_zscore_labels.py over {self.config.derived_root}"
                )
            frames.append(zscores)
        if not frames:
            return DataFrame(columns=["sample_id", "arm", "zscore"])
        return concat(frames, ignore_index=True)

    def assemble(self) -> DataFrame:
        samples = self.discover_samples()
        if not samples:
            raise ValueError(f"no completed samples found under {self.config.derived_root}")

        matrix = self.load_sample_metadata(samples)

        fragment_stats = self.load_long_table(samples, "fragment_stats_path")
        ichor_estimates = self.load_long_table(samples, "ichor_estimate_path")
        window_ratios = pivot_window_ratios(self.load_long_table(samples, "window_ratio_path")).reset_index()
        cnv_bins = pivot_cnv_bins(self.load_long_table(samples, "cnv_bin_path")).reset_index()
        arm_zscores = pivot_arm_zscores(self.load_arm_zscores_long(samples), column_key="arm").reset_index()

        # Left-joined onto the discovered sample list in §9's column order, so a sample
        # whose table is present but empty surfaces as an all-NaN row rather than
        # vanishing from the cohort.
        for other in (fragment_stats, arm_zscores, window_ratios, cnv_bins, ichor_estimates):
            matrix = matrix.merge(other, on="sample_id", how="left")
        return matrix

    def write(self, matrix: DataFrame) -> Path:
        self.config.output_path.parent.mkdir(parents=True, exist_ok=True)
        matrix.to_csv(self.config.output_path, index=False)
        return self.config.output_path

    def get_feature_columns(self, matrix: DataFrame) -> list[str]:
        """Feature columns -- everything that is not metadata or the ichorCNA estimate.

        The ichorCNA estimate is excluded for the same reason `feature_matrix` tags it
        ablation-only (phase0 §2): here it additionally doubles as the stand-in
        tumor-fraction label, so feeding it back in as a feature would be circular.
        """
        excluded = set(COHORT_METADATA_COLUMNS) | {COHORT_TUMOR_FRACTION_COLUMN}
        return [column for column in matrix.columns if column not in excluded]

    def summarise(self, matrix: DataFrame) -> dict[str, object]:
        """Cohort-level shape/completeness summary, for the run log and the QC report.

        `complete_feature_columns` is the number of feature columns present for *every*
        sample. It matters because per-sample window and arm sets differ, so the pivot
        is inherently ragged -- and the §6 decoy classifier can only use columns that
        survive that intersection.
        """
        feature_columns = self.get_feature_columns(matrix)
        features = matrix[feature_columns]
        return {
            "n_samples": int(len(matrix)),
            "n_feature_columns": len(feature_columns),
            "complete_feature_columns": int((features.notna().all(axis=0)).sum()),
            "samples_by_source_dataset": {
                str(dataset): int(count) for dataset, count in matrix["source_dataset"].value_counts().items()
            },
            "tumor_fraction_available": int(matrix[COHORT_TUMOR_FRACTION_COLUMN].notna().sum()),
            "tumor_fraction_nonzero": int((matrix[COHORT_TUMOR_FRACTION_COLUMN].fillna(0.0) > 0).sum()),
        }
