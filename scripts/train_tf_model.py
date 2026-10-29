"""Driver for the tumor-fraction regressor (phase0 §4).

Reads a Phase 3 §9 feature matrix CSV and writes one trained LightGBM model per
regime (`matched_pair`, `real_background` -- never pooled) under --output-root.
Every hyperparameter comes from workflow/config/tf_model.yaml via the eventual
Snakemake rule, not from defaults here, matching this repo's established
config-not-hardcoded pattern. Evaluation (RMSE, correlation vs. ichorCNA) is a
separate ticket/script; this one only trains and reports feature importances.
"""

from __future__ import annotations

import json
from argparse import ArgumentParser
from pathlib import Path

from pandas import read_csv

from flash_fraq.tf_model import FEATURE_SET_ALL, VALID_FEATURE_SETS, TFModelConfig, TumorFractionModelTrainer


def main() -> None:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--feature-matrix-path", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--importances-path", type=Path)
    parser.add_argument("--train-split-group", default="train")
    parser.add_argument("--validation-split-group", default="validation")
    parser.add_argument("--test-split-group", default="test")
    parser.add_argument(
        "--include-ichorcna-feature",
        action="store_true",
        help="Ablation variant that adds ichorCNA_tumor_fraction_estimate as a feature. "
        "Omit for the primary/headline model (phase0 §2).",
    )
    parser.add_argument(
        "--feature-set",
        choices=VALID_FEATURE_SETS,
        default=FEATURE_SET_ALL,
        help="'all' (default) trains on every §9 feature; 'cnv_only' restricts to the "
        "CNV/coverage block (cnv_logr_bin_*, arm_zscore_*, "
        "gc_corrected_total_coverage_window_*), the fallback when fragmentomics "
        "features are not interpretable in this cohort (audit finding 1/2); "
        "'fragment_only' restricts to the DELFI fragment-ratio block "
        "(short_long_ratio_window_*, fragment_size_*) -- cnv_only's complement, "
        "testing whether fragment-ratio signal alone recovers tumor fraction "
        "(docs/motivation-and-goal.md §1).",
    )
    parser.add_argument("--n-estimators", type=int, default=300)
    parser.add_argument("--num-leaves", type=int, default=31)
    parser.add_argument("--max-depth", type=int, default=-1)
    parser.add_argument("--learning-rate", type=float, default=0.05)
    parser.add_argument("--subsample", type=float, default=0.8)
    parser.add_argument("--colsample-bytree", type=float, default=0.8)
    parser.add_argument("--min-child-samples", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    config = TFModelConfig.from_sources(
        output_root=args.output_root,
        hyperparameters={
            "n_estimators": args.n_estimators,
            "num_leaves": args.num_leaves,
            "max_depth": args.max_depth,
            "learning_rate": args.learning_rate,
            "subsample": args.subsample,
            "colsample_bytree": args.colsample_bytree,
            "min_child_samples": args.min_child_samples,
            "random_state": args.seed,
        },
        train_split_group=args.train_split_group,
        validation_split_group=args.validation_split_group,
        test_split_group=args.test_split_group,
        include_ichorcna_feature=args.include_ichorcna_feature,
        feature_set=args.feature_set,
    )

    matrix = read_csv(args.feature_matrix_path)
    trainer = TumorFractionModelTrainer(config)
    models = trainer.run(matrix)

    report: dict[str, object] = {"variant": config.variant, "feature_set": config.feature_set, "regimes": {}}
    for regime, model in sorted(models.items()):
        if model is None:
            # audit finding 5: no row in this regime carries the training split
            # label -- no booster was fit, no joblib was written. A clear report
            # entry, not a silently missing key.
            report["regimes"][regime] = {
                "artifact_path": None,
                "n_training_rows": 0,
                "skipped_reason": f"no rows labeled split_group={args.train_split_group!r}; no booster fit",
            }
            print(f"{regime}: skipped -- no {args.train_split_group!r} rows to fit on")
            continue
        report["regimes"][regime] = {
            "artifact_path": str(config.model_artifact_path(regime)),
            "n_training_rows": model.n_training_rows,
            "top_features": model.feature_importances().head(10).to_dict("records"),
        }
        print(f"{regime}: trained on {model.n_training_rows} rows -> {config.model_artifact_path(regime)}")

    if args.importances_path is not None:
        args.importances_path.parent.mkdir(parents=True, exist_ok=True)
        args.importances_path.write_text(json.dumps(report, indent=2))

    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
