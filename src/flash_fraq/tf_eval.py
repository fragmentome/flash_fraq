"""Held-out evaluation + ichorCNA head-to-head comparison
(Phase 0 §5, docs/todo/todo/held-out-evaluation-and-ichorcna-comparison.md).

This is the project's "first number": it takes a `TrainedTumorFractionModel`
(`flash_fraq.tf_model`) that has already been fit, and *only* evaluates it -- it
never re-fits or re-tunes against held-out rows. Consumes the same §9-shaped
feature matrix the model was trained on, restricted to the `split_group` values
that were withheld from training (default: `validation` and `test`).

Per regime (`matched_pair`, `real_background` -- reported separately, never
pooled, mirroring `tf_model`'s own regime discipline) this module computes:

- **Primary metrics:** RMSE/MAE of the model's predictions vs.
  `empirical_tumor_fraction` (`REGRESSION_LABEL_COLUMN`).
- **ichorCNA baseline, side by side:** the *same* RMSE/MAE/correlation metrics
  computed for `ichorCNA_tumor_fraction_estimate` against the same label on the
  same rows, plus the correlation between the model's predictions and ichorCNA's
  estimate directly -- so "our model vs. ichorCNA" is a real comparison, not just
  an aggregate number for one side.
- **Detection floor:** the smallest true-tumor-fraction bin whose predictions are
  statistically distinguishable (one-sided Mann-Whitney U, `alternative="greater"`)
  from a null pool of `nominal_tumor_fraction == 0` rows, and remain so for every
  higher bin -- comparable to ichorCNA's own ~3% published floor (phase0 §1.1,
  §5.1). The null pool is defined on the mixture design's nominal fraction, not
  the empirical label: after the label estimator fix (audit finding 7), f=0
  mixtures still carry small positive empirical values from sequencing noise, so
  a label-based null pool (`empirical_tumor_fraction <= 0`) stays empty forever
  (audit finding 6).
- **Low-fraction sensitivity:** the same RMSE/MAE, restricted to rows below a
  configurable low-tumor-fraction threshold, for both the model and the ichorCNA
  baseline.
- **Disagreement analysis:** a row-level table of cases where the model and
  ichorCNA differ meaningfully in their error against the true label -- phase0
  §5.1's explicit "actual deliverable," not the aggregate correlation. HCC1395 is a
  fixed aneuploid line on both admixture regimes (`docs/admixture-design.md`), so
  neither regime currently has a near-diploid/low-CNA-burden axis; this module
  cannot characterize "ichorCNA-uninformative near-diploid" disagreement on data
  built this way, and surfaces that gap as a note rather than silently omitting it
  or pretending the check ran.
- **Feature importances**, passed through unchanged from
  `TrainedTumorFractionModel.feature_importances()`.

Split honesty: this module never re-derives or re-shuffles `split_group` --
rows are selected by their existing `split_group` label only (`eval_split_groups`,
default `("validation", "test")`), the same leakage-prevention key `tf_model`
already respects during `.fit()`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from pandas import DataFrame, Series
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import LeaveOneOut, cross_val_predict

from flash_fraq.feature_matrix import REGRESSION_LABEL_COLUMN
from flash_fraq.tf_model import FEATURE_SET_CNV_ONLY, REGIME_COLUMN, TrainedTumorFractionModel, get_feature_columns

ICHORCNA_COLUMN: str = "ichorCNA_tumor_fraction_estimate"
NOMINAL_FRACTION_COLUMN: str = "nominal_tumor_fraction"

# audit finding 5: matched_pair (11 rows, one per fraction rung, a single
# subject pair) has no train/validation/test split by construction -- see
# admixture_design.MATCHED_PAIR_SPLIT_GROUP and tf_model._train_rows, which now
# refuses to fit a booster for it at all. `compute_matched_pair_descriptive_report`
# is what this module reports instead: never a held-out evaluation.
MATCHED_PAIR_REGIME: str = "matched_pair"
MATCHED_PAIR_RIDGE_ALPHA: float = 1.0
MATCHED_PAIR_NOT_HELD_OUT_NOTE: str = (
    "matched_pair has one row per fraction rung, all from a single tumor/background "
    "subject pair, so it cannot be split by subject the way real_background is "
    "(admixture_design.MATCHED_PAIR_SPLIT_GROUP) and tf_model fits no booster for it "
    "(audit finding 5). Every number in this report is descriptive: label vs nominal "
    "fraction, ichorCNA vs label, and a leave-one-rung-out ridge fit on the CNV "
    "feature block. None of it is a held-out evaluation."
)

# Surfaced verbatim in every regime's report -- see the module docstring's
# "Disagreement analysis" note. Not a per-run computed fact (no near-diploid
# feature exists yet to detect it from); a fixed caveat about the current
# admixture design (docs/admixture-design.md), to be revisited once/if a
# near-diploid axis is added.
NO_NEAR_DIPLOID_AXIS_NOTE: str = (
    "Both admixture regimes (matched_pair, real_background) use HCC1395, a fixed "
    "aneuploid cell line, as the tumor source -- there is no near-diploid/low-CNA-"
    "burden axis in the current design. The 'ichorCNA is expected to be "
    "uninformative here' half of the disagreement analysis (phase0 §5.1) is "
    "therefore unsatisfiable as specified until such an axis exists; the "
    "disagreement table below characterizes model-vs-ichorCNA divergence in "
    "general, not specifically in the near-diploid regime this criterion asks for."
)


@dataclass(frozen=True)
class TFEvalConfig:
    """Everything `TumorFractionModelEvaluator` needs, independent of any one matrix.

    Every tunable here is declared in `workflow/config/tf_eval.yaml`, never
    hardcoded, matching this repo's established pattern (`tf_model.yaml`,
    `admixture_design.yaml`, ...).
    """

    eval_split_groups: tuple[str, ...] = ("validation", "test")
    low_fraction_threshold: float = 0.05
    disagreement_threshold: float = 0.1
    detection_floor_null_max: float = 0.0
    detection_floor_bin_edges: tuple[float, ...] = (
        0.01,
        0.02,
        0.03,
        0.05,
        0.08,
        0.12,
        0.2,
        0.3,
        0.4,
        0.5,
    )
    detection_floor_alpha: float = 0.05
    detection_floor_min_null_n: int = 5
    detection_floor_min_bin_n: int = 5

    @classmethod
    def from_sources(cls, **kwargs: object) -> "TFEvalConfig":
        defaults = cls()

        def _tuple_of_str(value: object, default: tuple[str, ...]) -> tuple[str, ...]:
            return tuple(str(v) for v in value) if value is not None else default

        def _tuple_of_float(value: object, default: tuple[float, ...]) -> tuple[float, ...]:
            return tuple(float(v) for v in value) if value is not None else default

        return cls(
            eval_split_groups=_tuple_of_str(kwargs.get("eval_split_groups"), defaults.eval_split_groups),
            low_fraction_threshold=float(kwargs.get("low_fraction_threshold", defaults.low_fraction_threshold)),
            disagreement_threshold=float(kwargs.get("disagreement_threshold", defaults.disagreement_threshold)),
            detection_floor_null_max=float(
                kwargs.get("detection_floor_null_max", defaults.detection_floor_null_max)
            ),
            detection_floor_bin_edges=_tuple_of_float(
                kwargs.get("detection_floor_bin_edges"), defaults.detection_floor_bin_edges
            ),
            detection_floor_alpha=float(kwargs.get("detection_floor_alpha", defaults.detection_floor_alpha)),
            detection_floor_min_null_n=int(
                kwargs.get("detection_floor_min_null_n", defaults.detection_floor_min_null_n)
            ),
            detection_floor_min_bin_n=int(
                kwargs.get("detection_floor_min_bin_n", defaults.detection_floor_min_bin_n)
            ),
        )


def get_eval_rows(matrix: DataFrame, regime: str, eval_split_groups: tuple[str, ...]) -> DataFrame:
    """Rows this evaluation is allowed to touch: one regime, restricted to the
    held-out `split_group` values -- never the training split, never re-derived."""
    if REGIME_COLUMN not in matrix.columns:
        raise ValueError(f"matrix has no {REGIME_COLUMN!r} column")
    regime_matrix = matrix[matrix[REGIME_COLUMN] == regime]
    if "split_group" not in regime_matrix.columns:
        raise ValueError(
            f"regime {regime!r} has no split_group column -- cannot select a held-out slice "
            "without re-deriving the split, which this module refuses to do"
        )
    return regime_matrix[regime_matrix["split_group"].isin(eval_split_groups)]


def _rmse(predictions: Series, truth: Series) -> float:
    return float(math.sqrt(((predictions.to_numpy() - truth.to_numpy()) ** 2).mean()))


def _mae(predictions: Series, truth: Series) -> float:
    return float((predictions.to_numpy() - truth.to_numpy()).__abs__().mean())


def _correlations(predictions: Series, truth: Series) -> dict[str, float | None]:
    """Pearson and Spearman correlation, both reported per the ticket's
    "Pearson/Spearman correlation" acceptance criterion -- not a single pick.
    `None` (rather than raising or NaN) when fewer than 2 rows or a constant
    series make correlation undefined."""
    if len(predictions) < 2 or predictions.nunique() < 2 or truth.nunique() < 2:
        return {"pearson": None, "spearman": None}
    pearson = float(pearsonr(predictions, truth)[0])
    spearman = float(spearmanr(predictions, truth)[0])
    return {"pearson": pearson, "spearman": spearman}


def compute_point_metrics(predictions: Series, truth: Series) -> dict[str, float | None]:
    """RMSE, MAE, and both correlations of `predictions` against `truth` -- the
    shared metric set used identically for the model and for the ichorCNA
    baseline, so the two are directly comparable."""
    metrics: dict[str, float | None] = {"rmse": _rmse(predictions, truth), "mae": _mae(predictions, truth)}
    metrics.update(_correlations(predictions, truth))
    return metrics


def compute_detection_floor(
    predictions: Series,
    truth: Series,
    config: TFEvalConfig,
    null_mask: Series | None = None,
) -> float | None:
    """Smallest true-tumor-fraction bin edge whose predictions are statistically
    distinguishable from the f≈0 null pool (Mann-Whitney U, one-sided:
    predictions greater than the null pool), and every higher bin also is.

    `null_mask` selects the null pool. When omitted this falls back to
    `truth <= config.detection_floor_null_max` -- but audit finding 6 found that,
    on the empirical label, no row at any rung has a truth value at or below 0
    even at nominal f=0 (small positive values from sequencing error against the
    alt allele), so that fallback pool stays permanently empty. Callers with a
    `nominal_tumor_fraction` column (`TumorFractionModelEvaluator.evaluate_regime`)
    pass `null_mask=nominal_tumor_fraction == 0` instead, which is immune to
    label noise since it is defined on the mixture design, not the label.

    Returns `None` when the null pool is too small to test against
    (`detection_floor_min_null_n`) or no bin (with enough rows,
    `detection_floor_min_bin_n`) clears `detection_floor_alpha` all the way up
    through the highest tested bin -- i.e. "undetectable given this data,"
    distinct from a floor of `0.0`.
    """
    if null_mask is None:
        null_mask = truth <= config.detection_floor_null_max
    null_predictions = predictions[null_mask]
    if len(null_predictions) < config.detection_floor_min_null_n:
        return None

    edges = list(config.detection_floor_bin_edges) + [math.inf]
    bin_results: list[tuple[float, float]] = []
    for lower, upper in zip(edges[:-1], edges[1:]):
        bin_mask = (truth >= lower) & (truth < upper)
        bin_predictions = predictions[bin_mask]
        if len(bin_predictions) < config.detection_floor_min_bin_n:
            continue
        _, p_value = mannwhitneyu(bin_predictions, null_predictions, alternative="greater")
        bin_results.append((lower, float(p_value)))

    for index, (lower, _) in enumerate(bin_results):
        if all(p < config.detection_floor_alpha for _, p in bin_results[index:]):
            return lower
    return None


def compute_disagreement_table(
    predictions: Series,
    ichorcna_estimate: Series,
    truth: Series,
    threshold: float,
    identity_columns: DataFrame | None = None,
) -> DataFrame:
    """Row-level model-vs-ichorCNA disagreement, not just an aggregate correlation
    (phase0 §5.1). A row is flagged when the two methods' absolute errors against
    the true label differ by more than `threshold`, labeled with which method was
    closer on that row."""
    model_abs_error = (predictions - truth).abs()
    ichorcna_abs_error = (ichorcna_estimate - truth).abs()
    error_gap = ichorcna_abs_error - model_abs_error  # positive: model closer

    table = DataFrame(
        {
            "true_tumor_fraction": truth,
            "model_prediction": predictions,
            "ichorcna_estimate": ichorcna_estimate,
            "model_abs_error": model_abs_error,
            "ichorcna_abs_error": ichorcna_abs_error,
            "error_gap_ichorcna_minus_model": error_gap,
            "model_closer": error_gap > 0,
        },
        index=predictions.index,
    )
    if identity_columns is not None:
        for column in identity_columns.columns:
            table.insert(0, column, identity_columns[column])

    disagreements = table[error_gap.abs() > threshold].sort_values(
        "error_gap_ichorcna_minus_model", ascending=False, key=lambda s: s.abs()
    )
    return disagreements.reset_index(drop=True)


def compute_leave_one_rung_out_ridge_predictions(
    features: DataFrame, label: Series, alpha: float = MATCHED_PAIR_RIDGE_ALPHA
) -> Series:
    """Leave-one-rung-out ridge predictions, sized for matched_pair's n=11 (one row
    per fraction rung): every prediction comes from a ridge fit on the other rungs,
    so this is a real -- if tiny -- generalisation check, never an in-sample fit.
    `LeaveOneOut` on 11 rows is exactly 11 one-row-out folds, one per rung.
    """
    predictions = cross_val_predict(Ridge(alpha=alpha, random_state=0), features, label, cv=LeaveOneOut())
    return Series(predictions, index=label.index, name="leave_one_rung_out_ridge_prediction")


@dataclass(frozen=True)
class MatchedPairDescriptiveReport:
    """Descriptive-only report for a regime too small to hold out by subject
    (matched_pair: audit finding 5). Never a held-out evaluation -- see
    `MATCHED_PAIR_NOT_HELD_OUT_NOTE`.
    """

    regime: str
    n_rows: int
    label_vs_nominal: dict[str, float | None]
    ichorcna_vs_label: dict[str, float | None]
    leave_one_rung_out_ridge_metrics: dict[str, float | None]
    n_cnv_features: int
    note: str = MATCHED_PAIR_NOT_HELD_OUT_NOTE

    def to_dict(self) -> dict[str, object]:
        return {
            "regime": self.regime,
            "n_rows": self.n_rows,
            "label_vs_nominal": self.label_vs_nominal,
            "ichorcna_vs_label": self.ichorcna_vs_label,
            "leave_one_rung_out_ridge_metrics": self.leave_one_rung_out_ridge_metrics,
            "n_cnv_features": self.n_cnv_features,
            "note": self.note,
        }


def compute_matched_pair_descriptive_report(
    matrix: DataFrame,
    regime: str = MATCHED_PAIR_REGIME,
    ridge_alpha: float = MATCHED_PAIR_RIDGE_ALPHA,
) -> MatchedPairDescriptiveReport:
    """Builds `regime`'s descriptive report directly from every row it has in
    `matrix` -- no split_group filtering, since matched_pair carries none (audit
    finding 5) and this report is explicitly not a held-out evaluation.
    """
    regime_matrix = matrix[matrix[REGIME_COLUMN] == regime]
    if regime_matrix.empty:
        raise ValueError(f"no rows for regime {regime!r}")
    if ICHORCNA_COLUMN not in regime_matrix.columns:
        raise ValueError(f"matrix is missing {ICHORCNA_COLUMN!r} -- required for the ichorCNA-vs-label metrics")
    if NOMINAL_FRACTION_COLUMN not in regime_matrix.columns:
        raise ValueError(f"matrix is missing {NOMINAL_FRACTION_COLUMN!r} -- required for the label-vs-nominal metrics")

    label = regime_matrix[REGRESSION_LABEL_COLUMN]
    nominal = regime_matrix[NOMINAL_FRACTION_COLUMN]
    ichorcna_estimate = regime_matrix[ICHORCNA_COLUMN]

    cnv_feature_columns = get_feature_columns(regime_matrix, feature_set=FEATURE_SET_CNV_ONLY)
    none_metrics = {"rmse": None, "mae": None, "pearson": None, "spearman": None}
    if cnv_feature_columns and len(regime_matrix) >= 3:
        ridge_predictions = compute_leave_one_rung_out_ridge_predictions(
            regime_matrix[cnv_feature_columns], label, alpha=ridge_alpha
        )
        ridge_metrics = compute_point_metrics(ridge_predictions, label)
    else:
        ridge_metrics = none_metrics

    return MatchedPairDescriptiveReport(
        regime=regime,
        n_rows=len(regime_matrix),
        label_vs_nominal=compute_point_metrics(label, nominal),
        ichorcna_vs_label=compute_point_metrics(ichorcna_estimate, label),
        leave_one_rung_out_ridge_metrics=ridge_metrics,
        n_cnv_features=len(cnv_feature_columns),
    )


@dataclass
class RegimeEvalResult:
    """One regime's full evaluation report -- never pooled with another regime's."""

    regime: str
    variant: str
    n_eval_rows: int
    model_metrics: dict[str, float | None]
    ichorcna_metrics: dict[str, float | None]
    model_vs_ichorcna_correlation: dict[str, float | None]
    low_fraction_threshold: float
    model_low_fraction_metrics: dict[str, float | None]
    ichorcna_low_fraction_metrics: dict[str, float | None]
    n_low_fraction_rows: int
    detection_floor: float | None
    disagreement: DataFrame
    feature_importances: DataFrame
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, object]:
        return {
            "regime": self.regime,
            "variant": self.variant,
            "n_eval_rows": self.n_eval_rows,
            "model_metrics": self.model_metrics,
            "ichorcna_metrics": self.ichorcna_metrics,
            "model_vs_ichorcna_correlation": self.model_vs_ichorcna_correlation,
            "low_fraction_threshold": self.low_fraction_threshold,
            "n_low_fraction_rows": self.n_low_fraction_rows,
            "model_low_fraction_metrics": self.model_low_fraction_metrics,
            "ichorcna_low_fraction_metrics": self.ichorcna_low_fraction_metrics,
            "detection_floor": self.detection_floor,
            "disagreement": self.disagreement.to_dict("records"),
            "feature_importances": self.feature_importances.to_dict("records"),
            "notes": self.notes,
        }


class TumorFractionModelEvaluator:
    """Evaluates already-trained `TrainedTumorFractionModel`s against their
    matrix's held-out split -- never fits, never tunes, never touches training
    rows."""

    def __init__(self, config: TFEvalConfig) -> None:
        self.config = config

    def evaluate_regime(self, model: TrainedTumorFractionModel, matrix: DataFrame) -> RegimeEvalResult:
        eval_rows = get_eval_rows(matrix, model.regime, self.config.eval_split_groups)
        if eval_rows.empty:
            raise ValueError(
                f"regime {model.regime!r} has no rows in split_group {self.config.eval_split_groups!r} "
                "to evaluate against"
            )
        if ICHORCNA_COLUMN not in eval_rows.columns:
            raise ValueError(f"matrix is missing {ICHORCNA_COLUMN!r} -- required for the ichorCNA comparison")
        if NOMINAL_FRACTION_COLUMN not in eval_rows.columns:
            raise ValueError(
                f"matrix is missing {NOMINAL_FRACTION_COLUMN!r} -- required for the detection-floor "
                "null pool (audit finding 6)"
            )

        truth = eval_rows[model.label_column]
        predictions = model.predict(eval_rows)
        ichorcna_estimate = eval_rows[ICHORCNA_COLUMN]
        # audit finding 6: the null pool is defined on the mixture design's nominal
        # fraction, not the empirical label -- see compute_detection_floor's docstring.
        detection_floor_null_mask = eval_rows[NOMINAL_FRACTION_COLUMN] == 0

        low_fraction_mask = truth < self.config.low_fraction_threshold
        low_fraction_rows = eval_rows[low_fraction_mask]
        low_fraction_predictions = predictions[low_fraction_mask]
        low_fraction_ichorcna = ichorcna_estimate[low_fraction_mask]
        low_fraction_truth = truth[low_fraction_mask]

        notes = [NO_NEAR_DIPLOID_AXIS_NOTE]

        return RegimeEvalResult(
            regime=model.regime,
            variant=model.variant,
            n_eval_rows=len(eval_rows),
            model_metrics=compute_point_metrics(predictions, truth),
            ichorcna_metrics=compute_point_metrics(ichorcna_estimate, truth),
            model_vs_ichorcna_correlation=_correlations(predictions, ichorcna_estimate),
            low_fraction_threshold=self.config.low_fraction_threshold,
            n_low_fraction_rows=len(low_fraction_rows),
            model_low_fraction_metrics=(
                compute_point_metrics(low_fraction_predictions, low_fraction_truth)
                if len(low_fraction_rows) >= 2
                else {"rmse": None, "mae": None, "pearson": None, "spearman": None}
            ),
            ichorcna_low_fraction_metrics=(
                compute_point_metrics(low_fraction_ichorcna, low_fraction_truth)
                if len(low_fraction_rows) >= 2
                else {"rmse": None, "mae": None, "pearson": None, "spearman": None}
            ),
            detection_floor=compute_detection_floor(
                predictions, truth, self.config, null_mask=detection_floor_null_mask
            ),
            disagreement=compute_disagreement_table(
                predictions,
                ichorcna_estimate,
                truth,
                self.config.disagreement_threshold,
                identity_columns=eval_rows[["sample_id"]] if "sample_id" in eval_rows.columns else None,
            ),
            feature_importances=model.feature_importances(),
            notes=notes,
        )

    def evaluate(
        self, models: dict[str, TrainedTumorFractionModel], matrix: DataFrame
    ) -> dict[str, RegimeEvalResult]:
        """One independent report per regime -- see class docstring on why regime
        results are never pooled or averaged together."""
        return {regime: self.evaluate_regime(model, matrix) for regime, model in models.items()}
