"""Driver for the held-out evaluation + ichorCNA comparison (phase0 §5).

Loads one trained `TrainedTumorFractionModel` per regime present in a §9 feature
matrix CSV (artifact paths matching `tf_model`'s `tf_model_{regime}_{variant}.joblib`
naming under --model-root) and evaluates each strictly against its held-out
`split_group` rows -- RMSE/MAE/correlation vs. the true label, the same metrics for
the ichorCNA baseline, model-vs-ichorCNA correlation, a low-fraction sensitivity
slice, a detection floor, and a row-level disagreement table. Training happens in
`scripts/train_tf_model.py`; this script only evaluates.
"""

from __future__ import annotations

import json
from argparse import ArgumentParser
from pathlib import Path

from pandas import read_csv

from flash_fraq.tf_eval import (
    MATCHED_PAIR_REGIME,
    TFEvalConfig,
    TumorFractionModelEvaluator,
    compute_matched_pair_descriptive_report,
)
from flash_fraq.tf_model import (
    FEATURE_SET_ALL,
    ICHORCNA_ABLATION_VARIANT,
    PRIMARY_MODEL_VARIANT,
    REGIME_COLUMN,
    VALID_FEATURE_SETS,
    TrainedTumorFractionModel,
)


def _model_artifact_path(model_root: Path, regime: str, variant: str, feature_set: str) -> Path:
    # Mirrors TFModelConfig.model_artifact_path's naming exactly -- kept as one
    # helper here rather than reconstructing a full TFModelConfig (which needs
    # hyperparameters this script never has) just to read that property.
    suffix = "" if feature_set == FEATURE_SET_ALL else f"_{feature_set}"
    return model_root / f"tf_model_{regime}_{variant}{suffix}.joblib"


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--feature-matrix-path", required=True, type=Path)
    parser.add_argument("--model-root", required=True, type=Path)
    parser.add_argument("--report-path", type=Path)
    parser.add_argument(
        "--variant",
        choices=(PRIMARY_MODEL_VARIANT, ICHORCNA_ABLATION_VARIANT),
        default=PRIMARY_MODEL_VARIANT,
        help="Which trained model variant to evaluate (default: the primary/headline model).",
    )
    parser.add_argument(
        "--feature-set",
        choices=VALID_FEATURE_SETS,
        default=FEATURE_SET_ALL,
        help="Which trained feature_set's model artifact to load (must match the "
        "--feature-set the model was trained with; see train_tf_model.py).",
    )
    parser.add_argument("--eval-split-groups", nargs="+", default=["validation", "test"])
    parser.add_argument("--low-fraction-threshold", type=float, default=0.05)
    parser.add_argument("--disagreement-threshold", type=float, default=0.1)
    parser.add_argument("--detection-floor-null-max", type=float, default=0.0)
    parser.add_argument(
        "--detection-floor-bin-edges",
        type=float,
        nargs="+",
        default=[0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.2, 0.3, 0.4, 0.5],
    )
    parser.add_argument("--detection-floor-alpha", type=float, default=0.05)
    parser.add_argument("--detection-floor-min-null-n", type=int, default=5)
    parser.add_argument("--detection-floor-min-bin-n", type=int, default=5)
    args = parser.parse_args()

    config = TFEvalConfig.from_sources(
        eval_split_groups=args.eval_split_groups,
        low_fraction_threshold=args.low_fraction_threshold,
        disagreement_threshold=args.disagreement_threshold,
        detection_floor_null_max=args.detection_floor_null_max,
        detection_floor_bin_edges=args.detection_floor_bin_edges,
        detection_floor_alpha=args.detection_floor_alpha,
        detection_floor_min_null_n=args.detection_floor_min_null_n,
        detection_floor_min_bin_n=args.detection_floor_min_bin_n,
    )

    matrix = read_csv(args.feature_matrix_path)
    if REGIME_COLUMN not in matrix.columns:
        raise ValueError(f"{args.feature_matrix_path}: matrix has no {REGIME_COLUMN!r} column")
    regimes = sorted(matrix[REGIME_COLUMN].unique())

    models = {}
    for regime in regimes:
        artifact_path = _model_artifact_path(args.model_root, regime, args.variant, args.feature_set)
        if not artifact_path.exists():
            # audit finding 5: matched_pair (and any other regime with no `train`
            # rows) never gets a joblib written by train_tf_model.py, so this is
            # the expected, non-error path for it -- the descriptive report below
            # is what covers it instead.
            print(f"skipping regime {regime!r}: no trained model at {artifact_path}")
            continue
        models[regime] = TrainedTumorFractionModel.load(artifact_path)

    if not models and MATCHED_PAIR_REGIME not in regimes:
        raise ValueError(f"no trained {args.variant} models found under {args.model_root} for regimes {regimes}")

    evaluator = TumorFractionModelEvaluator(config)
    results = evaluator.evaluate(models, matrix)

    report = {regime: result.to_dict() for regime, result in sorted(results.items())}
    for regime, result in sorted(results.items()):
        print(
            f"{regime}: n={result.n_eval_rows} "
            f"model rmse={result.model_metrics['rmse']:.4f} "
            f"ichorCNA rmse={result.ichorcna_metrics['rmse']:.4f} "
            f"detection_floor={result.detection_floor}"
        )

    if MATCHED_PAIR_REGIME in regimes:
        # audit finding 5: not held out by subject -- see
        # tf_eval.MATCHED_PAIR_NOT_HELD_OUT_NOTE, carried in the report itself.
        matched_pair_report = compute_matched_pair_descriptive_report(matrix)
        report[f"{MATCHED_PAIR_REGIME}_descriptive"] = matched_pair_report.to_dict()
        print(
            f"{MATCHED_PAIR_REGIME} (descriptive, not held-out): n={matched_pair_report.n_rows} "
            f"label_vs_nominal pearson={matched_pair_report.label_vs_nominal['pearson']} "
            f"ichorcna_vs_label pearson={matched_pair_report.ichorcna_vs_label['pearson']}"
        )

    if args.report_path is not None:
        args.report_path.parent.mkdir(parents=True, exist_ok=True)
        args.report_path.write_text(json.dumps(report, indent=2))

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
