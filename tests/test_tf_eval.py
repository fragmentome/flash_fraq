"""Held-out evaluation + ichorCNA comparison
(docs/todo/todo/held-out-evaluation-and-ichorcna-comparison.md).

Local-only per the ticket: this suite never trains against real data or a real
ichorCNA call. Two kinds of fixtures are used, deliberately:

- `synthetic_eval_matrix` extends `test_tf_model.py`'s §9-shaped fixture-building
  approach (same columns, same regime/split_group semantics) with an *injected*
  detection floor and a discrete near-zero "null" cluster, and trains a real
  `TrainedTumorFractionModel` via `flash_fraq.tf_model` on it, then evaluates
  that trained model -- this exercises the whole predict -> evaluate path the way
  a real run would.
- `test_compute_*` unit tests build tiny hand-computable `Series` directly,
  without any model, to verify each metric independently against a value worked
  out by hand -- "doesn't crash" is not enough per the ticket's own testing
  section, which asks for hand-computed expected values.

CORRELATION_FLOOR here matches test_tf_model.py's (0.6): same synthetic
signal-to-noise design, so the same bar applies to predictions on the held-out
slice.
"""

from pathlib import Path

import numpy
import pytest
from pandas import DataFrame, Series

from flash_fraq.tf_eval import (
    ICHORCNA_COLUMN,
    MATCHED_PAIR_NOT_HELD_OUT_NOTE,
    MATCHED_PAIR_REGIME,
    NO_NEAR_DIPLOID_AXIS_NOTE,
    TFEvalConfig,
    TumorFractionModelEvaluator,
    compute_detection_floor,
    compute_disagreement_table,
    compute_matched_pair_descriptive_report,
    compute_point_metrics,
    get_eval_rows,
)
from flash_fraq.tf_model import TFModelConfig, TumorFractionModelTrainer

REAL_BACKGROUND = "real_background"

INFORMATIVE_FEATURES = (
    "short_long_ratio_window_1",
    "gc_corrected_total_coverage_window_3",
    "arm_zscore_chr1p",
    "fragment_size_entropy",
)
N_NOISE_FEATURES = 24
N_WINDOWS = 10
ARM_LABELS = ("chr1p", "chr1q", "chr2p", "chr2q")
N_CNV_BINS = 6

# The detection-floor value deliberately injected into the synthetic held-out
# data below: predictions are indistinguishable from the f=0 null pool for true
# fractions below this, and cleanly distinguishable at and above it.
INJECTED_DETECTION_FLOOR = 0.08


def _base_columns(n_rows: int, rng: numpy.random.Generator, prefix: str) -> dict[str, object]:
    columns: dict[str, object] = {
        "sample_id": [f"{prefix}_{i:04d}" for i in range(n_rows)],
        "tumor_source_id": ["HCC1395_WGS_NV_3"] * n_rows,
        "platform_tumor": ["HiSeq4000"] * n_rows,
        "platform_healthy": ["HiSeq2500"] * n_rows,
        "reference_build": ["hg19"] * n_rows,
        "target_coverage": [1.0] * n_rows,
        "mtdna_fraction_log10": rng.normal(-2.2, 0.1, n_rows),
        "fragment_size_median": rng.normal(167.0, 1.0, n_rows),
        "fragment_size_mode": rng.normal(166.0, 1.0, n_rows),
        "fragment_size_entropy": rng.normal(0.0, 1.0, n_rows),
    }
    for window_index in range(1, N_WINDOWS + 1):
        columns[f"short_long_ratio_window_{window_index}"] = rng.normal(0.55, 0.02, n_rows)
        columns[f"gc_corrected_total_coverage_window_{window_index}"] = rng.normal(10.0, 0.3, n_rows)
        columns[f"gc_corrected_short_coverage_window_{window_index}"] = rng.normal(4.0, 0.2, n_rows)
    for arm in ARM_LABELS:
        columns[f"arm_zscore_{arm}"] = rng.normal(0.0, 1.0, n_rows)
    for bin_index in range(1, N_CNV_BINS + 1):
        columns[f"cnv_logr_bin_{bin_index}"] = rng.normal(0.0, 0.1, n_rows)
    for noise_index in range(N_NOISE_FEATURES):
        columns[f"noise_feature_{noise_index}"] = rng.normal(0.0, 1.0, n_rows)
    return columns


def _inject_signal(
    columns: dict[str, object],
    true_fraction: numpy.ndarray,
    rng: numpy.random.Generator,
    informative_features: tuple[str, ...],
    noise_sigma: float = 0.03,
) -> None:
    for weight, feature in enumerate(informative_features, start=1):
        columns[feature] = weight * true_fraction + rng.normal(0.0, noise_sigma, len(true_fraction))


def _split_labels(n_rows: int, rng: numpy.random.Generator) -> list[str]:
    return list(rng.choice(["train", "validation", "test"], size=n_rows, p=[0.7, 0.15, 0.15]))


def _make_matrix(n_rows: int, seed: int, prefix: str, true_fraction: numpy.ndarray) -> DataFrame:
    rng = numpy.random.default_rng(seed)
    columns = _base_columns(n_rows, rng, prefix)
    _inject_signal(columns, true_fraction, rng, INFORMATIVE_FEATURES)
    # A noisy, correlated-but-not-identical proxy for the label -- ichorCNA's
    # published behavior -- with more noise than the model's informative
    # features get, so the "model beats ichorCNA" comparison in the tests below
    # is a real, non-vacuous difference rather than an artifact of equal noise.
    columns[ICHORCNA_COLUMN] = true_fraction + rng.normal(0.0, 0.12, n_rows)
    columns["nominal_tumor_fraction"] = true_fraction
    columns["empirical_tumor_fraction"] = true_fraction
    columns["regime"] = REAL_BACKGROUND
    columns["healthy_source_id"] = [f"{prefix}_bg_{i:04d}" for i in range(n_rows)]
    columns["tumor_subject_id"] = ["HCC1395"] * n_rows
    columns["healthy_subject_id"] = [f"{prefix}_subject_{i:04d}" for i in range(n_rows)]
    columns["split_group"] = _split_labels(n_rows, rng)
    return DataFrame(columns)


@pytest.fixture
def synthetic_eval_matrix() -> DataFrame:
    """A single-regime (`real_background`) §9-shaped matrix, train/validation/test
    split, with true fractions concentrated below and above `INJECTED_DETECTION_FLOOR`
    (plus an explicit f=0 cluster) so the detection-floor routine has something
    real to recover."""
    rng = numpy.random.default_rng(99)
    n_rows = 500
    # A mix of: an explicit null cluster (f=0), a "below floor" band, and rows
    # spanning up to 0.5 -- rather than pure uniform -- so every detection-floor
    # bin this suite checks against has enough rows.
    zero_cluster = numpy.zeros(60)
    below_floor = rng.uniform(0.0, INJECTED_DETECTION_FLOOR, 120)
    above_floor = rng.uniform(INJECTED_DETECTION_FLOOR, 0.5, n_rows - 60 - 120)
    true_fraction = numpy.concatenate([zero_cluster, below_floor, above_floor])
    rng.shuffle(true_fraction)
    return _make_matrix(n_rows, seed=7, prefix="mix_rb", true_fraction=true_fraction)


def _train_config(output_root: Path) -> TFModelConfig:
    # tf_model.py's production defaults (300 estimators, shallow learning rate) --
    # used as-is here (rather than the smaller/faster hyperparameters
    # test_tf_model.py uses) because this fixture's tighter signal budget (four
    # informative features scaled 1x-4x true_fraction, noise_sigma=0.03) needs
    # more boosting rounds to cleanly separate informative features' importance
    # from the 24 pure-noise columns'.
    #
    # colsample_bytree is pinned to 1.0 here (production's default is 0.8): with
    # per-tree random feature subsampling, LightGBM's fixed random_state draws a
    # column subset keyed off the *total* feature count, so which noise column
    # gets over-selected is an accident of exactly how many feature columns
    # exist -- audit finding 2 correctly dropped mtdna_fraction_log10 from the
    # candidate feature set (get_feature_columns), which changed that count by
    # one and flipped this fixture's informative-vs-noise margin (verified: the
    # gap between the weakest informative feature and the strongest noise
    # feature is thin, and sign-flips across colsample_bytree's random subset
    # draw for several unrelated feature counts/seeds, not just this one). Using
    # every feature at every split removes that per-run subsampling accident
    # entirely, so this assertion tests the real signal-vs-noise separation
    # instead of an artifact of column count.
    return TFModelConfig.from_sources(output_root=output_root, hyperparameters={"colsample_bytree": 1.0})


def _eval_config(**overrides) -> TFEvalConfig:
    return TFEvalConfig.from_sources(**overrides)


# --- Unit tests: hand-computable metrics, no model involved -----------------


def test_compute_point_metrics_matches_hand_computed_values() -> None:
    predictions = Series([0.1, 0.2, 0.3, 0.4])
    truth = Series([0.0, 0.2, 0.5, 0.4])
    # errors: 0.1, 0.0, -0.2, 0.0 -> squared: 0.01, 0, 0.04, 0 -> mean 0.0125
    expected_rmse = 0.0125**0.5
    # abs errors: 0.1, 0.0, 0.2, 0.0 -> mean 0.075
    expected_mae = 0.075

    metrics = compute_point_metrics(predictions, truth)

    assert metrics["rmse"] == pytest.approx(expected_rmse)
    assert metrics["mae"] == pytest.approx(expected_mae)
    # Hand-computed Pearson r for these two vectors (verified independently):
    # cov(pred, truth) / (std(pred) * std(truth))
    pred_arr = predictions.to_numpy()
    truth_arr = truth.to_numpy()
    expected_pearson = numpy.corrcoef(pred_arr, truth_arr)[0, 1]
    assert metrics["pearson"] == pytest.approx(expected_pearson)


def test_compute_point_metrics_returns_none_correlation_for_constant_series() -> None:
    predictions = Series([0.2, 0.2, 0.2])
    truth = Series([0.1, 0.2, 0.3])

    metrics = compute_point_metrics(predictions, truth)

    assert metrics["pearson"] is None
    assert metrics["spearman"] is None
    assert metrics["rmse"] is not None


def test_disagreement_table_flags_only_rows_past_the_threshold_and_labels_who_is_closer() -> None:
    truth = Series([0.5, 0.5, 0.5], index=[0, 1, 2])
    predictions = Series([0.5, 0.5, 0.9], index=[0, 1, 2])  # perfect, perfect, way off
    ichorcna = Series([0.5, 0.9, 0.5], index=[0, 1, 2])  # perfect, way off, perfect

    table = compute_disagreement_table(predictions, ichorcna, truth, threshold=0.1)

    assert len(table) == 2
    row1 = table[table["model_abs_error"] == 0.0].iloc[0]
    assert row1["model_closer"]
    row2 = table[table["ichorcna_abs_error"] == 0.0].iloc[0]
    assert not row2["model_closer"]


def test_detection_floor_recovers_an_injected_floor_on_a_hand_built_series() -> None:
    rng = numpy.random.default_rng(0)
    # Null pool: true fraction 0, predictions centered at 0 with small noise.
    null_truth = numpy.zeros(30)
    null_predictions = rng.normal(0.0, 0.01, 30)
    # Below the injected floor (0.05): predictions still centered at ~0, i.e.
    # indistinguishable from the null pool by construction.
    below_truth = rng.uniform(0.0, 0.05, 60)
    below_predictions = rng.normal(0.0, 0.01, 60)
    # At/above the injected floor: predictions track the true fraction cleanly.
    above_truth = rng.uniform(0.05, 0.5, 200)
    above_predictions = above_truth + rng.normal(0.0, 0.01, 200)

    truth = Series(numpy.concatenate([null_truth, below_truth, above_truth]))
    predictions = Series(numpy.concatenate([null_predictions, below_predictions, above_predictions]))

    config = _eval_config(
        detection_floor_bin_edges=[0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2, 0.3, 0.4, 0.5]
    )
    floor = compute_detection_floor(predictions, truth, config)

    assert floor == pytest.approx(0.05)


def test_detection_floor_is_none_when_null_pool_too_small() -> None:
    truth = Series(numpy.linspace(0.0, 0.5, 20))
    predictions = Series(numpy.linspace(0.0, 0.5, 20))
    config = _eval_config(detection_floor_min_null_n=100)

    assert compute_detection_floor(predictions, truth, config) is None


def test_detection_floor_null_mask_overrides_the_truth_based_default() -> None:
    """audit finding 6: after the label fix, f=0 mixtures still carry small
    positive empirical values (sequencing noise against the alt allele), so
    `truth <= detection_floor_null_max` finds an empty null pool even though the
    mixture design's nominal fraction is exactly 0 for those rows."""
    rng = numpy.random.default_rng(1)
    n_null = 30
    n_above = 200
    null_truth = rng.uniform(1e-4, 1e-3, n_null)  # slightly positive, never <= 0
    null_predictions = rng.normal(0.0, 0.01, n_null)
    above_truth = rng.uniform(0.05, 0.5, n_above)
    above_predictions = above_truth + rng.normal(0.0, 0.01, n_above)

    truth = Series(numpy.concatenate([null_truth, above_truth]))
    predictions = Series(numpy.concatenate([null_predictions, above_predictions]))
    config = _eval_config(detection_floor_bin_edges=[0.05, 0.08, 0.12, 0.2, 0.3, 0.4, 0.5])

    # The truth-based default finds nothing at or below 0 -- the pool the old
    # `truth <= 0.0` null-pool definition (audit finding 6) stayed empty on.
    assert compute_detection_floor(predictions, truth, config) is None

    null_mask = Series([True] * n_null + [False] * n_above)
    assert compute_detection_floor(predictions, truth, config, null_mask=null_mask) is not None


# --- matched_pair descriptive report (not held-out) --------------------------


def _matched_pair_matrix(n_rungs: int = 11, seed: int = 11) -> DataFrame:
    """A matched_pair-shaped matrix: one row per fraction rung, CNV features that
    track nominal fraction, and a non-CNV column that must never reach the
    leave-one-rung-out ridge fit."""
    rng = numpy.random.default_rng(seed)
    nominal = numpy.linspace(0.0, 0.5, n_rungs)
    columns: dict[str, object] = {
        "sample_id": [f"mp_{i:02d}" for i in range(n_rungs)],
        "regime": MATCHED_PAIR_REGIME,
        "nominal_tumor_fraction": nominal,
        "empirical_tumor_fraction": nominal + rng.normal(0.0, 0.01, n_rungs),
        ICHORCNA_COLUMN: nominal + rng.normal(0.0, 0.05, n_rungs),
        "fragment_size_entropy": rng.normal(0.0, 1.0, n_rungs),  # non-CNV: excluded
    }
    for bin_index in range(1, 4):
        columns[f"cnv_logr_bin_{bin_index}"] = nominal * (0.1 * bin_index) + rng.normal(0.0, 0.01, n_rungs)
    for arm in ("chr1p", "chr1q"):
        columns[f"arm_zscore_{arm}"] = nominal * 2.0 + rng.normal(0.0, 0.05, n_rungs)
    return DataFrame(columns)


def test_matched_pair_descriptive_report_computes_label_vs_nominal_and_ichorcna_vs_label() -> None:
    matrix = _matched_pair_matrix()

    report = compute_matched_pair_descriptive_report(matrix)

    assert report.regime == MATCHED_PAIR_REGIME
    assert report.n_rows == 11
    assert report.label_vs_nominal["pearson"] is not None
    assert report.label_vs_nominal["pearson"] > 0.9
    assert report.ichorcna_vs_label["pearson"] is not None
    assert report.leave_one_rung_out_ridge_metrics["rmse"] is not None
    assert report.note == MATCHED_PAIR_NOT_HELD_OUT_NOTE


def test_matched_pair_descriptive_report_ridge_uses_only_cnv_features() -> None:
    matrix = _matched_pair_matrix()

    report = compute_matched_pair_descriptive_report(matrix)

    # 3 cnv_logr_bin_* + 2 arm_zscore_* -- fragment_size_entropy is non-CNV and
    # must never be counted or used.
    assert report.n_cnv_features == 5


def test_matched_pair_descriptive_report_raises_on_missing_nominal_fraction() -> None:
    matrix = _matched_pair_matrix().drop(columns=["nominal_tumor_fraction"])

    with pytest.raises(ValueError, match="nominal_tumor_fraction"):
        compute_matched_pair_descriptive_report(matrix)


# --- get_eval_rows: split honesty, no leakage --------------------------------


def test_get_eval_rows_only_returns_held_out_split_rows() -> None:
    matrix = DataFrame(
        {
            "regime": [REAL_BACKGROUND] * 6,
            "split_group": ["train", "train", "validation", "validation", "test", "test"],
            "value": [1, 2, 3, 4, 5, 6],
        }
    )

    eval_rows = get_eval_rows(matrix, REAL_BACKGROUND, ("validation", "test"))

    assert set(eval_rows["split_group"]) == {"validation", "test"}
    assert "train" not in set(eval_rows["split_group"])
    assert sorted(eval_rows["value"]) == [3, 4, 5, 6]


def test_get_eval_rows_respects_regime_partition() -> None:
    matrix = DataFrame(
        {
            "regime": ["real_background", "matched_pair"],
            "split_group": ["validation", "validation"],
            "value": [1, 2],
        }
    )

    eval_rows = get_eval_rows(matrix, "matched_pair", ("validation",))

    assert list(eval_rows["value"]) == [2]


# --- Full pipeline: train a real model, evaluate it, check no leakage -------


def test_evaluator_end_to_end_metrics_and_ichorcna_comparison_and_no_leakage(
    tmp_path: Path, synthetic_eval_matrix: DataFrame
) -> None:
    """Leakage canary: after training (on the clean matrix, so the model actually
    learns the real signal), build a corrupted *copy* of the matrix whose
    training-only rows' informative features and label are overwritten with a
    signal completely disjoint from the true-fraction relationship -- eval rows
    are left untouched. If `evaluate()` read anything from the training rows
    (pooled a correlation across all rows, computed a floor over the whole
    matrix, ...), the corruption would change its output. It does not: results
    on the corrupted matrix are byte-for-byte identical to results on the clean
    one, which is the actual no-leakage check, not just a row-count assertion."""
    matrix = synthetic_eval_matrix

    trainer = TumorFractionModelTrainer(_train_config(tmp_path))
    models = trainer.train(matrix)

    evaluator = TumorFractionModelEvaluator(_eval_config())
    clean_results = evaluator.evaluate(models, matrix)

    corrupted_matrix = matrix.copy()
    train_mask = corrupted_matrix["split_group"] == "train"
    for feature in INFORMATIVE_FEATURES:
        corrupted_matrix.loc[train_mask, feature] = 999.0
    corrupted_matrix.loc[train_mask, "empirical_tumor_fraction"] = 999.0
    corrupted_matrix.loc[train_mask, ICHORCNA_COLUMN] = 999.0

    corrupted_results = evaluator.evaluate(models, corrupted_matrix)

    result = clean_results[REAL_BACKGROUND]
    corrupted_result = corrupted_results[REAL_BACKGROUND]

    assert corrupted_result.model_metrics == result.model_metrics
    assert corrupted_result.ichorcna_metrics == result.ichorcna_metrics
    assert corrupted_result.detection_floor == result.detection_floor
    assert corrupted_result.n_eval_rows == result.n_eval_rows

    # Split honesty: the evaluator's own eval slice must be exactly the
    # validation+test rows, never train rows.
    expected_eval_rows = matrix[
        (matrix["regime"] == REAL_BACKGROUND) & (matrix["split_group"].isin(["validation", "test"]))
    ]
    assert result.n_eval_rows == len(expected_eval_rows)
    assert result.n_eval_rows < len(matrix)  # sanity: it's a strict subset

    # Primary metrics computed and sane on a signal this generous.
    assert result.model_metrics["rmse"] is not None
    assert result.model_metrics["pearson"] > 0.5

    # ichorCNA baseline computed independently, on the same rows, same label.
    assert result.ichorcna_metrics["rmse"] is not None
    assert result.ichorcna_metrics["pearson"] is not None

    # Hand-computable cross-check: ichorCNA's RMSE against truth on eval rows
    # equals what we'd get calling compute_point_metrics ourselves.
    eval_rows = get_eval_rows(matrix, REAL_BACKGROUND, ("validation", "test"))
    expected_ichorcna_metrics = compute_point_metrics(
        eval_rows[ICHORCNA_COLUMN], eval_rows["empirical_tumor_fraction"]
    )
    assert result.ichorcna_metrics["rmse"] == pytest.approx(expected_ichorcna_metrics["rmse"])
    assert result.ichorcna_metrics["mae"] == pytest.approx(expected_ichorcna_metrics["mae"])

    # ichorCNA's estimate here was built noisier than the model's informative
    # features (0.12 sigma vs 0.03), so the model should out-RMSE it on this
    # fixture -- a real, checkable "our model vs ichorCNA" comparison.
    assert result.model_metrics["rmse"] < result.ichorcna_metrics["rmse"]

    # Detection floor recovered close to the injected floor.
    assert result.detection_floor is not None
    assert result.detection_floor <= 0.12  # comfortably within the injected-floor neighborhood

    # Feature importances passed through, real feature names, every informative
    # feature ranked above every noise feature.
    importances = result.feature_importances.set_index("feature")["importance"]
    noise_columns = [c for c in matrix.columns if c.startswith("noise_feature_")]
    assert importances.loc[list(INFORMATIVE_FEATURES)].min() > importances.loc[noise_columns].max()

    # Disagreement table is a real DataFrame with the documented columns.
    for column in ("true_tumor_fraction", "model_prediction", "ichorcna_estimate", "model_closer"):
        assert column in result.disagreement.columns

    # The near-diploid-axis gap is surfaced, not silently dropped.
    assert NO_NEAR_DIPLOID_AXIS_NOTE in result.notes

    # Low-fraction sensitivity slice computed on the configured threshold.
    assert result.n_low_fraction_rows == int(
        (expected_eval_rows["empirical_tumor_fraction"] < result.low_fraction_threshold).sum()
    )


def test_evaluator_raises_on_missing_ichorcna_column(tmp_path: Path, synthetic_eval_matrix: DataFrame) -> None:
    trainer = TumorFractionModelTrainer(_train_config(tmp_path))
    models = trainer.train(synthetic_eval_matrix)

    stripped = synthetic_eval_matrix.drop(columns=[ICHORCNA_COLUMN])
    evaluator = TumorFractionModelEvaluator(_eval_config())

    with pytest.raises(ValueError, match="ichorCNA_tumor_fraction_estimate"):
        evaluator.evaluate(models, stripped)


def test_evaluator_raises_when_eval_split_has_no_rows(tmp_path: Path, synthetic_eval_matrix: DataFrame) -> None:
    trainer = TumorFractionModelTrainer(_train_config(tmp_path))
    models = trainer.train(synthetic_eval_matrix)

    evaluator = TumorFractionModelEvaluator(_eval_config(eval_split_groups=["nonexistent_split"]))

    with pytest.raises(ValueError, match="nonexistent_split"):
        evaluator.evaluate(models, synthetic_eval_matrix)
