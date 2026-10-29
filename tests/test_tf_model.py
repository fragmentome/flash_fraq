"""Tumor-fraction regressor (docs/todo/todo/tumor-fraction-regressor.md).

Local-only per the ticket: a synthetic §9-shaped feature matrix
(flash_fraq.feature_matrix's schema -- window ratios, arm z-scores, CNV bins,
fragment-size stats, mtdna fraction, ichorCNA estimate, empirical_tumor_fraction
label, split_group/regime) with a known, injectable tumor-fraction signal planted
in a handful of named features and pure noise everywhere else. The real §9 matrix
does not exist yet (feature-matrix-assembler is on-hold); this suite is the ticket's
entire close bar.

CORRELATION_FLOOR = 0.6 (Pearson r on the held-out `validation` split, real_background
regime): chosen as a clear, comfortably-clearable bar for a signal-to-noise ratio this
generous (four informative features, one dominant) trained on a few hundred synthetic
rows -- not a claim about real-data performance, which is the next ticket's job.
"""

from pathlib import Path

import numpy
import pytest
from pandas import DataFrame

from flash_fraq.tf_model import (
    ICHORCNA_ABLATION_VARIANT,
    PRIMARY_MODEL_VARIANT,
    TFModelConfig,
    TrainedTumorFractionModel,
    TumorFractionModelTrainer,
    get_feature_columns,
)

REAL_BACKGROUND = "real_background"
MATCHED_PAIR = "matched_pair"

# Feature names shaped like the §9 schema (flash_fraq.feature_matrix /
# docs/todo/on-hold/feature-matrix-assembler.md), not generic column names --
# the ticket requires importances reported with real feature names.
INFORMATIVE_FEATURES_REAL_BACKGROUND = (
    "short_long_ratio_window_1",
    "gc_corrected_total_coverage_window_3",
    "arm_zscore_chr1p",
    "fragment_size_entropy",
)
INFORMATIVE_FEATURES_MATCHED_PAIR = (
    "short_long_ratio_window_7",
    "arm_zscore_chr2q",
)
N_NOISE_FEATURES = 24
N_WINDOWS = 10
ARM_LABELS = ("chr1p", "chr1q", "chr2p", "chr2q")
N_CNV_BINS = 6


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
        # ablation-only per phase0 §2 -- deliberately made *informative* of the
        # label here (a plausible estimate would correlate with true tumor
        # fraction) so the exclusion tests are meaningful rather than vacuous.
        "ichorCNA_tumor_fraction_estimate": None,
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
    noise_sigma: float = 0.05,
) -> None:
    """Overwrite `informative_features` with true_fraction * weight + small noise;
    every other numeric column (already written by `_base_columns`) stays pure noise,
    independent of `true_fraction` by construction."""
    n_rows = len(true_fraction)
    for weight, feature in enumerate(informative_features, start=1):
        columns[feature] = weight * true_fraction + rng.normal(0.0, noise_sigma, n_rows)


def _make_regime_matrix(
    n_rows: int,
    seed: int,
    prefix: str,
    regime: str,
    informative_features: tuple[str, ...],
    split_group: object,
) -> DataFrame:
    rng = numpy.random.default_rng(seed)
    true_fraction = rng.uniform(0.0, 0.5, n_rows)
    columns = _base_columns(n_rows, rng, prefix)
    _inject_signal(columns, true_fraction, rng, informative_features)
    # A noisy but real proxy for the label, the way an actual ichorCNA point
    # estimate would be -- correlated, not identical, and only allowed to help
    # when a caller opts into the ablation variant.
    columns["ichorCNA_tumor_fraction_estimate"] = true_fraction + rng.normal(0.0, 0.08, n_rows)
    columns["nominal_tumor_fraction"] = true_fraction
    columns["empirical_tumor_fraction"] = true_fraction
    columns["regime"] = regime
    columns["healthy_source_id"] = [f"{prefix}_bg_{i:04d}" for i in range(n_rows)]
    columns["tumor_subject_id"] = ["HCC1395"] * n_rows
    columns["healthy_subject_id"] = [f"{prefix}_subject_{i:04d}" for i in range(n_rows)]
    if isinstance(split_group, str):
        columns["split_group"] = [split_group] * n_rows
    else:
        columns["split_group"] = split_group
    return DataFrame(columns)


def _split_labels(n_rows: int, rng: numpy.random.Generator) -> list[str]:
    # 70/15/15 train/validation/test, matching admixture_design.yaml's split_fractions.
    return list(rng.choice(["train", "validation", "test"], size=n_rows, p=[0.7, 0.15, 0.15]))


@pytest.fixture
def synthetic_matrix() -> DataFrame:
    """A two-regime §9-shaped matrix: `real_background` is split (train/validation/
    test, like the real design); `matched_pair` is unsplit (split_group == "",
    mirroring admixture_design.MATCHED_PAIR_SPLIT_GROUP) and carries a *different*
    informative feature set, so a test can tell the two regimes' models apart."""
    n_real_background = 400
    rng = numpy.random.default_rng(1234)
    real_background = _make_regime_matrix(
        n_real_background,
        seed=1,
        prefix="mix_rb",
        regime=REAL_BACKGROUND,
        informative_features=INFORMATIVE_FEATURES_REAL_BACKGROUND,
        split_group=_split_labels(n_real_background, rng),
    )
    matched_pair = _make_regime_matrix(
        60,
        seed=2,
        prefix="mix_mp",
        regime=MATCHED_PAIR,
        informative_features=INFORMATIVE_FEATURES_MATCHED_PAIR,
        split_group="",  # unsplit, as the real design declares matched_pair
    )
    matrix = DataFrame(
        list(real_background.to_dict("records")) + list(matched_pair.to_dict("records"))
    )
    return matrix


def _config(output_root: Path, **overrides) -> TFModelConfig:
    hyperparameters = {
        "n_estimators": 80,
        "num_leaves": 7,
        "max_depth": -1,
        "learning_rate": 0.15,
        "subsample": 0.9,
        "colsample_bytree": 0.9,
        "min_child_samples": 5,
        "random_state": 42,
    }
    hyperparameters.update(overrides.pop("hyperparameters", {}))
    return TFModelConfig.from_sources(output_root=output_root, hyperparameters=hyperparameters, **overrides)


CORRELATION_FLOOR = 0.6


def test_feature_importances_rank_informative_features_above_noise(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    trainer = TumorFractionModelTrainer(_config(tmp_path))
    models = trainer.train(synthetic_matrix)

    importances = models[REAL_BACKGROUND].feature_importances()
    top_features = set(importances.head(len(INFORMATIVE_FEATURES_REAL_BACKGROUND))["feature"])

    assert top_features == set(INFORMATIVE_FEATURES_REAL_BACKGROUND)

    informative_importance = importances.set_index("feature").loc[list(INFORMATIVE_FEATURES_REAL_BACKGROUND), "importance"]
    noise_columns = [c for c in synthetic_matrix.columns if c.startswith("noise_feature_")]
    noise_importance = importances.set_index("feature").loc[noise_columns, "importance"]
    assert informative_importance.min() > noise_importance.max()


def test_predictions_correlate_with_true_fraction_on_held_out_split(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    trainer = TumorFractionModelTrainer(_config(tmp_path))
    models = trainer.train(synthetic_matrix)

    regime_matrix = synthetic_matrix[synthetic_matrix["regime"] == REAL_BACKGROUND]
    held_out = regime_matrix[regime_matrix["split_group"] == "validation"]
    assert len(held_out) > 10  # sanity: the fixture actually produced a held-out slice

    predictions = models[REAL_BACKGROUND].predict(held_out)
    correlation = numpy.corrcoef(predictions, held_out["empirical_tumor_fraction"])[0, 1]

    assert correlation > CORRELATION_FLOOR


def test_regimes_are_trained_separately_not_pooled(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    """real_background's informative features must be found on their own -- if
    matched_pair's disjoint informative set (which never reaches `.fit()`, see
    `test_matched_pair_regime_with_no_train_rows_is_skipped`) leaked in through
    pooling, it would show up in this ranking too."""
    trainer = TumorFractionModelTrainer(_config(tmp_path))
    models = trainer.train(synthetic_matrix)

    rb_top = set(models[REAL_BACKGROUND].feature_importances().head(len(INFORMATIVE_FEATURES_REAL_BACKGROUND))["feature"])

    assert rb_top == set(INFORMATIVE_FEATURES_REAL_BACKGROUND)
    # Neither regime's informative set leaked into the other's top ranking, which
    # is what silently pooling the two regimes would produce.
    assert rb_top.isdisjoint(INFORMATIVE_FEATURES_MATCHED_PAIR)


def test_matched_pair_regime_with_no_train_rows_is_skipped(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    """matched_pair carries split_group == "" on every row (admixture_design's
    MATCHED_PAIR_SPLIT_GROUP sentinel) -- since no row is ever labeled "train",
    the trainer must not fit a booster for it at all (audit finding 5: the old
    fallback to training on every row produced an 11-row constant predictor)."""
    trainer = TumorFractionModelTrainer(_config(tmp_path))
    models = trainer.train(synthetic_matrix)

    assert set(models) == {REAL_BACKGROUND, MATCHED_PAIR}
    assert models[REAL_BACKGROUND] is not None
    assert models[MATCHED_PAIR] is None


def test_save_all_writes_no_artifact_for_a_skipped_regime(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    trainer = TumorFractionModelTrainer(_config(tmp_path))
    models = trainer.run(synthetic_matrix)

    assert models[MATCHED_PAIR] is None
    assert trainer.config.model_artifact_path(REAL_BACKGROUND).exists()
    assert not trainer.config.model_artifact_path(MATCHED_PAIR).exists()


def test_real_background_training_only_uses_the_train_split(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    """The trainer must consume split_group as written, never re-shuffle rows: the
    number of rows actually fit on must equal exactly the rows labeled 'train'."""
    trainer = TumorFractionModelTrainer(_config(tmp_path))
    models = trainer.train(synthetic_matrix)

    regime_matrix = synthetic_matrix[synthetic_matrix["regime"] == REAL_BACKGROUND]
    n_train_rows = int((regime_matrix["split_group"] == "train").sum())
    assert models[REAL_BACKGROUND].n_training_rows == n_train_rows
    assert 0 < n_train_rows < len(regime_matrix)


def test_split_group_is_respected_no_subject_leakage_across_split(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    """Every source subject in this fixture is unique to one row (like the real
    admixture design's background pool at the subject level), so this checks the
    mechanism directly: predicting on a row from one split must never have used
    that row's own label during training, which holds iff the trained model's
    row count matches the train-only split exactly (checked above) and the
    val/test subjects are disjoint from the train subjects by construction here."""
    regime_matrix = synthetic_matrix[synthetic_matrix["regime"] == REAL_BACKGROUND]
    train_subjects = set(regime_matrix.loc[regime_matrix["split_group"] == "train", "healthy_subject_id"])
    other_subjects = set(regime_matrix.loc[regime_matrix["split_group"] != "train", "healthy_subject_id"])
    assert train_subjects.isdisjoint(other_subjects)


def test_ichorcna_feature_excluded_by_default(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    trainer = TumorFractionModelTrainer(_config(tmp_path))
    models = trainer.train(synthetic_matrix)

    assert models[REAL_BACKGROUND].variant == PRIMARY_MODEL_VARIANT
    assert "ichorCNA_tumor_fraction_estimate" not in models[REAL_BACKGROUND].feature_columns


def test_ichorcna_ablation_variant_includes_the_feature_when_requested(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    config = _config(tmp_path, include_ichorcna_feature=True)
    trainer = TumorFractionModelTrainer(config)
    models = trainer.train(synthetic_matrix)

    assert models[REAL_BACKGROUND].variant == ICHORCNA_ABLATION_VARIANT
    assert "ichorCNA_tumor_fraction_estimate" in models[REAL_BACKGROUND].feature_columns
    # The two variants must be distinctly labeled artifact paths -- never the same
    # file -- so the ablation can never silently overwrite the headline model.
    assert config.model_artifact_path(REAL_BACKGROUND) != _config(tmp_path).model_artifact_path(REAL_BACKGROUND)


def test_get_feature_columns_excludes_identity_label_and_partition_columns(synthetic_matrix: DataFrame) -> None:
    feature_columns = get_feature_columns(synthetic_matrix, include_ichorcna_feature=False)

    for excluded in (
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
        "empirical_tumor_fraction",
        "split_group",
        "regime",
        "ichorCNA_tumor_fraction_estimate",
        "mtdna_fraction_log10",
    ):
        assert excluded not in feature_columns
    assert "short_long_ratio_window_1" in feature_columns


def test_get_feature_columns_cnv_only_selects_only_the_cnv_block(synthetic_matrix: DataFrame) -> None:
    feature_columns = get_feature_columns(synthetic_matrix, feature_set="cnv_only")

    assert feature_columns  # sanity: the fixture actually has CNV/coverage columns
    for column in feature_columns:
        assert column.startswith(("cnv_logr_bin_", "arm_zscore_", "gc_corrected_total_coverage_window_"))
    # Excluded blocks: short/long ratios, short coverage, fragment-size stats,
    # mtDNA fraction -- none of these are in CNV_ONLY_FEATURE_PREFIXES.
    for excluded_example in (
        "short_long_ratio_window_1",
        "gc_corrected_short_coverage_window_1",
        "fragment_size_entropy",
        "mtdna_fraction_log10",
        "noise_feature_0",
    ):
        assert excluded_example not in feature_columns


def test_get_feature_columns_fragment_only_selects_only_the_fragment_ratio_block(
    synthetic_matrix: DataFrame,
) -> None:
    feature_columns = get_feature_columns(synthetic_matrix, feature_set="fragment_only")

    assert feature_columns  # sanity: the fixture actually has fragment-ratio columns
    for column in feature_columns:
        assert column.startswith(("short_long_ratio_window_", "fragment_size_"))
    # Excluded blocks: CNV bins, arm z-scores, total *and* short coverage (short
    # coverage still tracks copy number -- see FRAGMENT_ONLY_FEATURE_PREFIXES'S
    # comment), mtDNA fraction -- none of these are in FRAGMENT_ONLY_FEATURE_PREFIXES.
    for excluded_example in (
        "cnv_logr_bin_1",
        "arm_zscore_chr1p",
        "gc_corrected_total_coverage_window_1",
        "gc_corrected_short_coverage_window_1",
        "mtdna_fraction_log10",
        "noise_feature_0",
    ):
        assert excluded_example not in feature_columns


def test_get_feature_columns_rejects_an_unknown_feature_set(synthetic_matrix: DataFrame) -> None:
    with pytest.raises(ValueError, match="feature_set"):
        get_feature_columns(synthetic_matrix, feature_set="bogus")


def test_trained_model_never_carries_mtdna_fraction_log10_as_a_feature(
    tmp_path: Path, synthetic_matrix: DataFrame
) -> None:
    """audit finding 2: mtdna_fraction_log10 is a tumor-read counter (HCC1395 gDNA
    ~10^-2.6 vs cfDNA ~10^-4.9), not a fragmentomics feature -- it must never appear
    in a saved model's feature_columns."""
    trainer = TumorFractionModelTrainer(_config(tmp_path))
    trainer.run(synthetic_matrix)

    loaded = TrainedTumorFractionModel.load(trainer.config.model_artifact_path(REAL_BACKGROUND))
    assert "mtdna_fraction_log10" not in loaded.feature_columns


def test_cnv_only_feature_set_trains_and_is_recorded_on_the_model(
    tmp_path: Path, synthetic_matrix: DataFrame
) -> None:
    config = _config(tmp_path, feature_set="cnv_only")
    trainer = TumorFractionModelTrainer(config)
    models = trainer.train(synthetic_matrix)

    model = models[REAL_BACKGROUND]
    assert model.feature_set == "cnv_only"
    for column in model.feature_columns:
        assert column.startswith(("cnv_logr_bin_", "arm_zscore_", "gc_corrected_total_coverage_window_"))
    # cnv_only gets a distinct artifact path from the default "all" model, so it
    # can never silently overwrite the primary model for the same regime/variant.
    assert config.model_artifact_path(REAL_BACKGROUND) != _config(tmp_path).model_artifact_path(REAL_BACKGROUND)


def test_fragment_only_feature_set_trains_and_is_recorded_on_the_model(
    tmp_path: Path, synthetic_matrix: DataFrame
) -> None:
    config = _config(tmp_path, feature_set="fragment_only")
    trainer = TumorFractionModelTrainer(config)
    models = trainer.train(synthetic_matrix)

    model = models[REAL_BACKGROUND]
    assert model.feature_set == "fragment_only"
    for column in model.feature_columns:
        assert column.startswith(("short_long_ratio_window_", "fragment_size_"))
    # fragment_only gets a distinct artifact path from both the default "all"
    # model and cnv_only, so none of the three can ever silently overwrite
    # another for the same regime/variant.
    assert config.model_artifact_path(REAL_BACKGROUND) != _config(tmp_path).model_artifact_path(REAL_BACKGROUND)
    assert config.model_artifact_path(REAL_BACKGROUND) != _config(
        tmp_path, feature_set="cnv_only"
    ).model_artifact_path(REAL_BACKGROUND)


def test_feature_set_rejected_at_config_construction(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="feature_set"):
        _config(tmp_path, feature_set="bogus")


def test_train_raises_without_a_regime_column() -> None:
    matrix = DataFrame({"empirical_tumor_fraction": [0.1, 0.2], "feature_1": [1.0, 2.0]})
    trainer = TumorFractionModelTrainer(_config(Path("/tmp/unused-tf-model-test")))

    with pytest.raises(ValueError, match="regime"):
        trainer.regimes(matrix)


def test_save_and_load_round_trip(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    trainer = TumorFractionModelTrainer(_config(tmp_path))
    models = trainer.run(synthetic_matrix)

    artifact_path = trainer.config.model_artifact_path(REAL_BACKGROUND)
    assert artifact_path.exists()
    # Never written into the repo -- it lands under the caller-given tmp_path.
    assert tmp_path in artifact_path.parents

    loaded = TrainedTumorFractionModel.load(artifact_path)
    regime_matrix = synthetic_matrix[synthetic_matrix["regime"] == REAL_BACKGROUND]
    held_out = regime_matrix[regime_matrix["split_group"] == "validation"]

    original_predictions = models[REAL_BACKGROUND].predict(held_out)
    reloaded_predictions = loaded.predict(held_out)
    assert numpy.allclose(original_predictions.to_numpy(), reloaded_predictions.to_numpy())


def test_predict_raises_on_missing_feature_columns(tmp_path: Path, synthetic_matrix: DataFrame) -> None:
    trainer = TumorFractionModelTrainer(_config(tmp_path))
    models = trainer.train(synthetic_matrix)

    incomplete = synthetic_matrix.drop(columns=["short_long_ratio_window_1"])
    with pytest.raises(ValueError, match="short_long_ratio_window_1"):
        models[REAL_BACKGROUND].predict(incomplete)
