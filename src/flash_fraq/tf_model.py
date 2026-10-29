"""Tumor-fraction regressor (Phase 0 §4, docs/todo/todo/tumor-fraction-regressor.md).

Trains a gradient-boosting regression model (LightGBM -- see the module-level note
below for why over XGBoost) mapping the Phase 3 §9 feature matrix
(`flash_fraq.feature_matrix`) to `empirical_tumor_fraction`. This module owns
*training* only: RMSE/correlation-vs-ichorCNA reporting and the source-held-out
generalisation check are a different, downstream ticket
(held-out-evaluation-and-ichorcna-comparison.md) that consumes the public interface
here (`TrainedTumorFractionModel.load` / `.predict` / `.feature_importances`) --
keep that interface stable.

Two design decisions this module exists to enforce, both already settled upstream
and restated here because getting either wrong would silently invalidate the model:

- **Regime is never pooled.** `docs/admixture-design.md`'s "Reopened decision"
  declares `matched_pair` and `real_background` two regimes evaluated separately,
  never merged into one training table. This module honours that by training one
  independent booster per regime (`TumorFractionModelTrainer.train`) rather than by
  passing `regime` in as an ordinary covariate -- a covariate a boosted tree is free
  to ignore if the loss doesn't need it, which is exactly the silent-pooling failure
  mode the admixture design forbids.
- **`split_group` is consumed, not re-derived.** Each regime's booster is fit only
  on rows whose `split_group` equals the configured training split; validation/test
  rows are never touched during `.fit()`. `matched_pair` carries no split (one
  background subject makes subject-level splitting impossible by construction, see
  `admixture_design.MATCHED_PAIR_SPLIT_GROUP`), so a regime whose rows never carry
  the training split label has nothing safe to fit on: `TumorFractionModelTrainer`
  skips it entirely (no booster, no joblib artifact) rather than falling back to
  training on every row -- audit finding 5, which found the old fallback produced
  an 11-row constant predictor with every feature importance at 0. `train()`'s
  result for such a regime is `None`; `TumorFractionModelTrainer._train_rows`
  documents when this triggers.

ichorCNA's point estimate (phase0 §2: "including it directly muddies the
complementary-to-ichorCNA story if it dominates") is excluded from the feature set
by default (`TFModelConfig.include_ichorcna_feature = False`, the "primary" model
variant) and only ever added back through an explicit, clearly-labeled ablation
variant (`include_ichorcna_feature=True`, `TrainedTumorFractionModel.variant ==
"ichorcna_ablation"`) -- never silently in the headline model.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

import joblib
from lightgbm import LGBMRegressor
from pandas import DataFrame, Series

from flash_fraq.feature_matrix import ABLATION_ONLY_COLUMNS, REGRESSION_LABEL_COLUMN

# LightGBM over XGBoost: both satisfy phase0 §4 (gradient boosting, not deep
# learning); LightGBM's PyPI wheel is a ~3.5 MB pure-CPU package, while
# xgboost>=2 pulls in an NVIDIA CUDA/NCCL wheel (~250 MB) as a runtime dependency
# even for CPU-only training on this small-N tabular problem -- pure bloat here.
# Both expose sklearn-compatible `.fit`/`.predict`/`.feature_importances_`, so
# nothing downstream is XGBoost/LightGBM-specific.

# Columns that describe a row's identity/provenance/label rather than a
# fragmentomics measurement -- never candidate features, regardless of dtype.
# Mirrors `feature_matrix.FeatureMatrixAssembler.get_primary_feature_columns`'s
# `non_feature_columns`, plus `regime` (a hard partition key this module trains
# separately on, never a feature within a regime's own model) and `split_group`
# (the leakage-prevention key, likewise never a predictor).
NON_FEATURE_COLUMNS: frozenset[str] = frozenset(
    {
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
        # admixture.AdmixtureRunner.write_sample_metadata sets REGRESSION_LABEL_COLUMN to
        # the *first* cross-estimator label's own empirical_tumor_fraction_<label> value
        # verbatim, then writes every label's value under its own column too -- so
        # empirical_tumor_fraction_sites (the canonical label in this cohort) is not
        # merely correlated with the regression label, it *is* it, and
        # empirical_tumor_fraction_clonal is the alternate estimator on the same SNVs.
        # Both are cross-estimator agreement diagnostics, never candidate features.
        "empirical_tumor_fraction_sites",
        "empirical_tumor_fraction_clonal",
        "split_group",
        "regime",
        # audit finding 2: HCC1395 gDNA carries ~10^-2.6 mtDNA fraction vs cfDNA's
        # ~10^-4.9, so this is a tumor-read counter (rho=0.96 vs nominal f in
        # real_background, #2 feature by split importance) rather than a
        # fragmentomics measurement -- stays in the matrix as metadata only.
        "mtdna_fraction_log10",
    }
)

REGIME_COLUMN: str = "regime"

PRIMARY_MODEL_VARIANT: str = "primary"
ICHORCNA_ABLATION_VARIANT: str = "ichorcna_ablation"

# Feature-set options (audit finding 2). "all" is every non-excluded numeric §9
# column, as before; "cnv_only" restricts to the CNV/coverage block -- the
# fallback for a cohort where the fragmentomics features are not interpretable
# (audit finding 1: sheared gDNA vs cfDNA library-prep confound). Mirrors
# `decoy_qc.CNV_FEATURE_PREFIXES` -- kept as its own tuple here rather than a
# shared import so this module's feature selection never depends on the QC
# gate's config.
FEATURE_SET_ALL: str = "all"
FEATURE_SET_CNV_ONLY: str = "cnv_only"
FEATURE_SET_FRAGMENT_ONLY: str = "fragment_only"
VALID_FEATURE_SETS: tuple[str, ...] = (FEATURE_SET_ALL, FEATURE_SET_CNV_ONLY, FEATURE_SET_FRAGMENT_ONLY)
CNV_ONLY_FEATURE_PREFIXES: tuple[str, ...] = (
    "cnv_logr_bin_",
    "arm_zscore_",
    "gc_corrected_total_coverage_window_",
)
# "fragment_only" (docs/motivation-and-goal.md §1: does DELFI-style fragment-ratio
# signal alone recover tumor fraction?) is cnv_only's complement -- the
# fragmentomics features with no CNV/coverage signal in them. Restricted to
# `short_long_ratio_window_*` (the DELFI fragment-length ratio itself) and
# `fragment_size_*` (median/mode/entropy of the fragment-length distribution).
# Deliberately excludes `gc_corrected_short_coverage_window_*`: it is a
# *coverage* quantity (short-fragment read depth per window), and per
# `CNV_ONLY_FEATURE_PREFIXES` above, GC-corrected coverage still tracks copy
# number -- including it here would leak CNV signal back into what is supposed
# to be the pure-fragmentomics feature set. Also excludes `mtdna_fraction_log10`
# (never a feature at all -- see `NON_FEATURE_COLUMNS`'s audit finding 2 note, a
# tumor-read counter, not a fragmentomics measurement) and
# `ichorcna_tumor_fraction_estimate` (phase0 §2's ablation-only feature, gated by
# `include_ichorcna_feature` same as every other feature set).
FRAGMENT_ONLY_FEATURE_PREFIXES: tuple[str, ...] = (
    "short_long_ratio_window_",
    "fragment_size_",
)


def get_feature_columns(
    matrix: DataFrame,
    include_ichorcna_feature: bool = False,
    feature_set: str = FEATURE_SET_ALL,
) -> list[str]:
    """Numeric §9 feature columns, real names, in schema column order.

    Excludes identity/label/partition columns (`NON_FEATURE_COLUMNS`) and, unless
    `include_ichorcna_feature` is set, `feature_matrix.ABLATION_ONLY_COLUMNS` --
    phase0 §2's ablation-only feature. Restricting to numeric dtypes is a cheap
    guard against a stray metadata string column being fed to the booster as a
    "feature" by accident.

    `feature_set="cnv_only"` (audit finding 2) further restricts the result to
    `CNV_ONLY_FEATURE_PREFIXES` -- the CNV bins, arm z-scores and GC-corrected
    total-coverage windows -- excluding every fragmentomics/short-fraction
    feature whose signal is confounded by tumor-side library prep in the
    current admixture design (audit finding 1).

    `feature_set="fragment_only"` (docs/motivation-and-goal.md §1) is cnv_only's
    complement: `FRAGMENT_ONLY_FEATURE_PREFIXES` -- the DELFI short/long ratio
    windows and fragment-size summary stats -- excluding every CNV/coverage
    feature (see that tuple's comment for why GC-corrected short coverage is
    excluded too).
    """
    if feature_set not in VALID_FEATURE_SETS:
        raise ValueError(f"feature_set must be one of {VALID_FEATURE_SETS}, got {feature_set!r}")

    excluded = set(NON_FEATURE_COLUMNS)
    if not include_ichorcna_feature:
        excluded |= set(ABLATION_ONLY_COLUMNS)
    columns = [
        column
        for column in matrix.columns
        if column not in excluded and _is_numeric(matrix[column])
    ]
    if feature_set == FEATURE_SET_CNV_ONLY:
        columns = [column for column in columns if column.startswith(CNV_ONLY_FEATURE_PREFIXES)]
    elif feature_set == FEATURE_SET_FRAGMENT_ONLY:
        columns = [column for column in columns if column.startswith(FRAGMENT_ONLY_FEATURE_PREFIXES)]
    return columns


def _is_numeric(series: Series) -> bool:
    return series.dtype.kind in "biuf"


@dataclass(frozen=True)
class TFModelHyperparameters:
    """LightGBM hyperparameters, declared in `workflow/config/tf_model.yaml` --
    never hardcoded here, per this repo's established pattern for every other
    tunable (the fraction ladder, `C_target`, ... in `admixture_design.yaml`)."""

    n_estimators: int = 300
    num_leaves: int = 31
    max_depth: int = -1
    learning_rate: float = 0.05
    subsample: float = 0.8
    colsample_bytree: float = 0.8
    min_child_samples: int = 10
    random_state: int = 42

    @classmethod
    def from_sources(cls, **kwargs: object) -> "TFModelHyperparameters":
        defaults = cls()
        return cls(
            n_estimators=int(kwargs.get("n_estimators", defaults.n_estimators)),
            num_leaves=int(kwargs.get("num_leaves", defaults.num_leaves)),
            max_depth=int(kwargs.get("max_depth", defaults.max_depth)),
            learning_rate=float(kwargs.get("learning_rate", defaults.learning_rate)),
            subsample=float(kwargs.get("subsample", defaults.subsample)),
            colsample_bytree=float(kwargs.get("colsample_bytree", defaults.colsample_bytree)),
            min_child_samples=int(kwargs.get("min_child_samples", defaults.min_child_samples)),
            random_state=int(kwargs.get("random_state", defaults.random_state)),
        )

    def as_lgbm_kwargs(self) -> dict[str, object]:
        return {
            "n_estimators": self.n_estimators,
            "num_leaves": self.num_leaves,
            "max_depth": self.max_depth,
            "learning_rate": self.learning_rate,
            "subsample": self.subsample,
            "colsample_bytree": self.colsample_bytree,
            "min_child_samples": self.min_child_samples,
            "random_state": self.random_state,
            # Silences LightGBM's per-tree stdout ("No further splits with
            # positive gain") which is routine noise on small/synthetic data
            # and not a training failure.
            "verbosity": -1,
        }


@dataclass(frozen=True)
class TFModelConfig:
    """Everything `TumorFractionModelTrainer` needs, independent of any one matrix."""

    hyperparameters: TFModelHyperparameters
    # Derived-data output root for serialized model artifacts -- weights are not
    # source, so this must never default into the repo tree. Callers (the CLI,
    # tests) always pass an explicit path.
    output_root: Path
    train_split_group: str = "train"
    validation_split_group: str = "validation"
    test_split_group: str = "test"
    include_ichorcna_feature: bool = False
    # audit finding 2's feature_set option -- "all" (default), "cnv_only", or
    # "fragment_only". Validated eagerly in __post_init__ so a typo surfaces at
    # config-construction time, not three hours into a training run.
    feature_set: str = FEATURE_SET_ALL

    def __post_init__(self) -> None:
        if self.feature_set not in VALID_FEATURE_SETS:
            raise ValueError(f"feature_set must be one of {VALID_FEATURE_SETS}, got {self.feature_set!r}")

    @property
    def variant(self) -> str:
        return ICHORCNA_ABLATION_VARIANT if self.include_ichorcna_feature else PRIMARY_MODEL_VARIANT

    @classmethod
    def from_sources(
        cls,
        output_root: Path,
        hyperparameters: dict[str, object] | None = None,
        train_split_group: str = "train",
        validation_split_group: str = "validation",
        test_split_group: str = "test",
        include_ichorcna_feature: bool = False,
        feature_set: str = FEATURE_SET_ALL,
    ) -> "TFModelConfig":
        return cls(
            hyperparameters=TFModelHyperparameters.from_sources(**(hyperparameters or {})),
            output_root=Path(output_root),
            train_split_group=train_split_group,
            validation_split_group=validation_split_group,
            test_split_group=test_split_group,
            include_ichorcna_feature=bool(include_ichorcna_feature),
            feature_set=str(feature_set),
        )

    def model_artifact_path(self, regime: str) -> Path:
        # The default feature_set keeps the pre-existing filename (never breaks an
        # existing artifact path); cnv_only and fragment_only each get their own
        # distinct suffix (derived from the feature_set name itself, so a new
        # feature_set needs no change here) so neither can ever silently
        # overwrite the primary/all-features model for the same regime and
        # variant.
        suffix = "" if self.feature_set == FEATURE_SET_ALL else f"_{self.feature_set}"
        return self.output_root / f"tf_model_{regime}_{self.variant}{suffix}.joblib"


@dataclass
class TrainedTumorFractionModel:
    """One fitted booster for one regime/variant, plus what a caller needs to use it.

    Public interface for the downstream evaluation ticket:

    - `predict(matrix)` -- tumor-fraction predictions for any §9-shaped matrix that
      carries this model's `feature_columns` (validation/test/external rows all
      work; nothing about `predict` depends on `split_group` or `regime`).
    - `feature_importances()` -- ranked `DataFrame` with real feature names (phase0
      §4's "built-in feature importance for the complementary-signal story").
    - `save(path)` / `TrainedTumorFractionModel.load(path)` -- round-trip a
      artifact through an explicit, caller-given derived-data path.
    """

    regime: str
    variant: str
    booster: LGBMRegressor
    feature_columns: tuple[str, ...]
    hyperparameters: TFModelHyperparameters
    train_split_group: str
    n_training_rows: int
    label_column: str = REGRESSION_LABEL_COLUMN
    feature_set: str = FEATURE_SET_ALL

    def predict(self, matrix: DataFrame) -> Series:
        """Tumor-fraction predictions, indexed like `matrix`."""
        missing = [column for column in self.feature_columns if column not in matrix.columns]
        if missing:
            raise ValueError(f"matrix is missing feature columns this model was trained on: {missing}")
        predictions = self.booster.predict(matrix[list(self.feature_columns)])
        return Series(predictions, index=matrix.index, name="predicted_tumor_fraction")

    def feature_importances(self) -> DataFrame:
        """Ranked feature importances (descending), real §9 feature names.

        LightGBM's default importance type ("split": how often a feature is used
        to split, summed over all trees) is what `LGBMRegressor` reports via
        `feature_importances_`; ties broken by feature name for determinism.
        """
        importances = DataFrame(
            {"feature": self.feature_columns, "importance": self.booster.feature_importances_}
        )
        return importances.sort_values(
            ["importance", "feature"], ascending=[False, True], ignore_index=True
        )

    def save(self, path: Path) -> Path:
        """Serialize this model under `path` -- never a repo path (weights are not source)."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(self, path)
        return path

    @staticmethod
    def load(path: Path) -> "TrainedTumorFractionModel":
        model = joblib.load(Path(path))
        if not isinstance(model, TrainedTumorFractionModel):
            raise TypeError(f"{path}: expected a TrainedTumorFractionModel, got {type(model)!r}")
        return model


class TumorFractionModelTrainer:
    """Fits one `TrainedTumorFractionModel` per regime present in a §9 matrix."""

    def __init__(self, config: TFModelConfig) -> None:
        self.config = config

    def regimes(self, matrix: DataFrame) -> list[str]:
        if REGIME_COLUMN not in matrix.columns:
            raise ValueError(
                f"matrix has no {REGIME_COLUMN!r} column -- the two-regime design "
                "(docs/admixture-design.md) requires it on every row"
            )
        return sorted(matrix[REGIME_COLUMN].unique())

    def _train_rows(self, regime_matrix: DataFrame) -> DataFrame | None:
        """Rows this regime's booster is allowed to fit on, or `None` when the
        regime has no row labeled the training split at all.

        Consumes `split_group` as written -- never re-shuffles rows -- so no
        source subject placed in validation/test by upstream split assignment can
        ever reach `.fit()`. A regime with no `split_group` column, or none of
        whose rows carry the training split label (e.g. `matched_pair`, unsplit
        by construction -- see `admixture_design.MATCHED_PAIR_SPLIT_GROUP`) has
        nothing safe to fit on: it used to fall back to training on every row,
        which on `matched_pair`'s 11 rows produced a booster that could not split
        under `min_child_samples=10` -- a constant predictor with every feature
        importance at 0 (audit finding 5). Returning `None` here tells
        `train_regime` to skip fitting entirely rather than reporting that
        constant predictor as a real model.
        """
        if "split_group" not in regime_matrix.columns:
            return None
        train_rows = regime_matrix[regime_matrix["split_group"] == self.config.train_split_group]
        if train_rows.empty:
            return None
        return train_rows

    def train_regime(self, matrix: DataFrame, regime: str) -> TrainedTumorFractionModel | None:
        """Fits and returns `regime`'s model, or `None` when it has no `train` rows
        to fit on (see `_train_rows`) -- in which case no booster is created and
        no artifact is ever written for this regime by `save_all`/`run`."""
        regime_matrix = matrix[matrix[REGIME_COLUMN] == regime]
        if regime_matrix.empty:
            raise ValueError(f"no rows for regime {regime!r}")

        train_rows = self._train_rows(regime_matrix)
        if train_rows is None:
            return None

        feature_columns = get_feature_columns(
            matrix, self.config.include_ichorcna_feature, self.config.feature_set
        )
        if not feature_columns:
            raise ValueError("no numeric feature columns found after excluding identity/label/partition columns")

        booster = LGBMRegressor(**self.config.hyperparameters.as_lgbm_kwargs())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            booster.fit(train_rows[feature_columns], train_rows[REGRESSION_LABEL_COLUMN])

        return TrainedTumorFractionModel(
            regime=regime,
            variant=self.config.variant,
            booster=booster,
            feature_columns=tuple(feature_columns),
            hyperparameters=self.config.hyperparameters,
            train_split_group=self.config.train_split_group,
            n_training_rows=len(train_rows),
            feature_set=self.config.feature_set,
        )

    def train(self, matrix: DataFrame) -> dict[str, TrainedTumorFractionModel | None]:
        """One independent model per regime -- see the module docstring on why
        regime is a training-loop partition, never a covariate handed to a single
        pooled booster. A regime with no `train` rows maps to `None` (audit
        finding 5) rather than being silently dropped from the result or fit on
        a fallback of every row."""
        return {regime: self.train_regime(matrix, regime) for regime in self.regimes(matrix)}

    def save_all(self, models: dict[str, TrainedTumorFractionModel | None]) -> dict[str, Path]:
        """Writes an artifact for every regime that actually has a model --
        `None` entries (no `train` rows, see `train_regime`) are skipped, so a
        regime that could not be fit never gets a joblib on disk."""
        return {
            regime: model.save(self.config.model_artifact_path(regime))
            for regime, model in models.items()
            if model is not None
        }

    def run(self, matrix: DataFrame) -> dict[str, TrainedTumorFractionModel | None]:
        """Train every regime's model and write its artifact under `config.output_root`."""
        models = self.train(matrix)
        self.save_all(models)
        return models
