"""Reshaped-mixture downstream rerun (docs/todo/todo/reshaped-mixture-downstream-rerun.md):
feature-matrix assembly -> decoy QC (per regime) -> tf_model (per feature_set) ->
tf_eval (per feature_set), over the flat mix_<name>/ admixture-mixture layout.

Deliberately NOT part of workflow/Snakefile / workflow/rules/local_features.smk:
see workflow/config/admixture_downstream_reshaped.yaml's header for why the
four per-mixture stages (relabel, fragments, GC-correct, ichorCNA) are run by
the existing bash drivers instead, outside this DAG entirely. This file only
models the four whole-matrix stages, each already a single tested script
(scripts/run_admixture_mixture_feature_matrix.py, scripts/run_decoy_qc.py,
scripts/train_tf_model.py, scripts/run_tf_eval.py) -- Snakemake's value here is
resumability (skip a stage whose output is newer than its input) and a real
`snakemake -n` job count per stage, not fan-out over hundreds of samples.

Included from workflow/Snakefile_reshaped_downstream, never from the main
workflow/Snakefile -- it must never share a DAG with the live
snakemake_workdirs/admixture build.
"""

from pathlib import Path

# Explicit interpreter, not a bare "python": scripts/drivers/reshaped_downstream_chain.sh
# runs snakemake itself from the "snakemake" conda env (same pattern as
# admixture_mixtures_local_driver.sh), which has snakemake's own deps but not
# flash_fraq's (pandas/lightgbm/scikit-learn) -- confirmed live: `conda activate
# snakemake && python -c "import sklearn"` raises ModuleNotFoundError. Every
# shell: command below that needs those packages (feature-matrix assembly,
# decoy QC, tf_model, tf_eval -- unlike local_features.smk's mosdepth-only
# rules, this file trains/evaluates with lightgbm+scikit-learn) must therefore
# invoke the flash_fraq env's interpreter by full path rather than relying on
# PATH/conda activate at snakemake-invocation time. The path itself lives in
# workflow/config/admixture_downstream_reshaped.yaml (not hardcoded here) so
# `snakemake --lint` doesn't flag an absolute path baked into the workflow.
FLASHFRAQ_PYTHON = config["flashfraq_python"]

DOWNSTREAM_CONFIG = config["admixture_downstream"]
DOWNSTREAM_MIXTURES_ROOT = Path(DOWNSTREAM_CONFIG["mixtures_root"])
DOWNSTREAM_COMBINED_PLAN = Path(DOWNSTREAM_CONFIG["combined_plan_path"])
DOWNSTREAM_SPLIT_ASSIGNMENT = Path(DOWNSTREAM_CONFIG["split_assignment_path"])
DOWNSTREAM_AGG_DIR = Path(DOWNSTREAM_CONFIG["agg_dir"])
DOWNSTREAM_FEATURE_MATRIX = Path(DOWNSTREAM_CONFIG["feature_matrix_path"])
DOWNSTREAM_FEATURE_MATRIX_SUMMARY = Path(DOWNSTREAM_CONFIG["feature_matrix_summary_path"])
DOWNSTREAM_FEATURE_MATRIX_REAL_BACKGROUND = Path(DOWNSTREAM_CONFIG["feature_matrix_real_background_path"])
DOWNSTREAM_MASK_MARGIN_BP = DOWNSTREAM_CONFIG.get("mask_margin_bp", 5000000)

DOWNSTREAM_DECOY_QC_DIR = Path(DOWNSTREAM_CONFIG["decoy_qc_dir"])
DOWNSTREAM_DECOY_QC_REGIMES = DOWNSTREAM_CONFIG["decoy_qc_regimes"]

DOWNSTREAM_MODELS_ROOT = Path(DOWNSTREAM_CONFIG["models_root"])
DOWNSTREAM_FEATURE_SETS = DOWNSTREAM_CONFIG["feature_sets"]
DOWNSTREAM_HYPERPARAMS = DOWNSTREAM_CONFIG["tf_model_hyperparameters"]


wildcard_constraints:
    regime="|".join(DOWNSTREAM_DECOY_QC_REGIMES),
    feature_set="|".join(DOWNSTREAM_FEATURE_SETS),


rule downstream_feature_matrix:
    """scripts/run_admixture_mixture_feature_matrix.py over every mix_* dir under
    DOWNSTREAM_MIXTURES_ROOT, joined against DOWNSTREAM_COMBINED_PLAN for `regime`
    (see that file's own comment for why it must be the *combined*, not the bare
    reshaped, plan). Also writes the real_background-only subset for parity with
    the original 495-mixture run's feature_matrix_real_background.csv.
    """
    input:
        combined_plan=str(DOWNSTREAM_COMBINED_PLAN),
        split_assignment=str(DOWNSTREAM_SPLIT_ASSIGNMENT),
    output:
        feature_matrix=str(DOWNSTREAM_FEATURE_MATRIX),
        summary=str(DOWNSTREAM_FEATURE_MATRIX_SUMMARY),
        real_background=str(DOWNSTREAM_FEATURE_MATRIX_REAL_BACKGROUND),
    params:
        mixtures_root=lambda wc, output: str(Path(output.feature_matrix).parent),
        agg_dir=str(DOWNSTREAM_AGG_DIR),
        mask_margin_bp=DOWNSTREAM_MASK_MARGIN_BP,
    log:
        str(DOWNSTREAM_MIXTURES_ROOT / "logs" / "downstream_feature_matrix.log"),
    conda:
        "flash_fraq"
    threads: 4
    resources:
        mem_mb=16384,
        runtime=120,
    shell:
        f"{FLASHFRAQ_PYTHON} {SCRIPTS_DIR / 'run_admixture_mixture_feature_matrix.py'} "
        "--mixtures-root {params.mixtures_root} "
        "--split-assignment-path {input.split_assignment} "
        "--admixture-plan-path {input.combined_plan} "
        "--agg-dir {params.agg_dir} "
        "--output-path {output.feature_matrix} "
        "--summary-path {output.summary} "
        "--mask-margin-bp {params.mask_margin_bp} "
        "> {log} 2>&1 && "
        f"{FLASHFRAQ_PYTHON} -c \""
        "import pandas as pd, sys; "
        "m = pd.read_csv(sys.argv[1]); "
        "m[m['regime'] == 'real_background'].to_csv(sys.argv[2], index=False)"
        "\" {output.feature_matrix} {output.real_background} >> {log} 2>&1"


rule downstream_decoy_qc:
    """Per-regime decoy-classifier gate (phase2 §6). Always run with
    --no-fail-on-flag so Snakemake's own output-on-failure deletion never
    removes a written report -- the HARD gate on `passed` is enforced by
    scripts/drivers/reshaped_downstream_chain.sh reading this report's JSON
    (decoy_qc.py's own report_path()/`passed` semantics), not by this rule's
    exit code.
    """
    input:
        feature_matrix=str(DOWNSTREAM_FEATURE_MATRIX),
    output:
        report=str(DOWNSTREAM_DECOY_QC_DIR / "decoy_qc_report_{regime}.json"),
    params:
        output_dir=lambda wc, output: str(Path(output.report).parent),
    log:
        str(DOWNSTREAM_MIXTURES_ROOT / "logs" / "downstream_decoy_qc" / "{regime}.log"),
    conda:
        "flash_fraq"
    threads: 2
    resources:
        mem_mb=8192,
        runtime=60,
    shell:
        f"{FLASHFRAQ_PYTHON} {SCRIPTS_DIR / 'run_decoy_qc.py'} "
        "--feature-matrix-path {input.feature_matrix} "
        "--output-dir {params.output_dir} "
        "--matrix-kind admixture "
        "--regime {wildcards.regime} "
        "--no-fail-on-flag "
        "> {log} 2>&1"


rule downstream_tf_model:
    """scripts/train_tf_model.py, primary variant (no ichorCNA ablation feature),
    one booster per regime present in the matrix -- never pooled. Gated on decoy
    QC passing by the orchestrator script, not by a Snakemake input edge (the
    report's `passed` field is a go/no-go the orchestrator reads, not a file this
    rule can depend on without coupling a hard business decision into the DAG).
    """
    input:
        feature_matrix=str(DOWNSTREAM_FEATURE_MATRIX),
    output:
        importances=str(DOWNSTREAM_MODELS_ROOT / "feature_importances_{feature_set}.json"),
    params:
        output_root=lambda wc, output: str(Path(output.importances).parent),
        n_estimators=DOWNSTREAM_HYPERPARAMS["n_estimators"],
        num_leaves=DOWNSTREAM_HYPERPARAMS["num_leaves"],
        max_depth=DOWNSTREAM_HYPERPARAMS["max_depth"],
        learning_rate=DOWNSTREAM_HYPERPARAMS["learning_rate"],
        subsample=DOWNSTREAM_HYPERPARAMS["subsample"],
        colsample_bytree=DOWNSTREAM_HYPERPARAMS["colsample_bytree"],
        min_child_samples=DOWNSTREAM_HYPERPARAMS["min_child_samples"],
        seed=DOWNSTREAM_HYPERPARAMS["random_state"],
    log:
        str(DOWNSTREAM_MIXTURES_ROOT.parent / "models" / "tf_model_reshaped" / "logs" / "{feature_set}.log"),
    conda:
        "flash_fraq"
    threads: 4
    resources:
        mem_mb=8192,
        runtime=60,
    shell:
        f"mkdir -p $(dirname {{log}}) && {FLASHFRAQ_PYTHON} {SCRIPTS_DIR / 'train_tf_model.py'} "
        "--feature-matrix-path {input.feature_matrix} "
        "--output-root {params.output_root} "
        "--importances-path {output.importances} "
        "--feature-set {wildcards.feature_set} "
        "--n-estimators {params.n_estimators} "
        "--num-leaves {params.num_leaves} "
        "--max-depth {params.max_depth} "
        "--learning-rate {params.learning_rate} "
        "--subsample {params.subsample} "
        "--colsample-bytree {params.colsample_bytree} "
        "--min-child-samples {params.min_child_samples} "
        "--seed {params.seed} "
        "> {log} 2>&1"


rule downstream_tf_eval:
    """scripts/run_tf_eval.py, primary variant, against the same feature matrix.
    Depends on this feature_set's feature_importances.json as a stand-in for
    "the model(s) for this feature_set have been trained" (train_tf_model.py
    writes both the regime joblib(s) and the importances report in the same
    invocation; matched_pair has no joblib by design -- see cnv-only-interim-model.md --
    but the importances report is still written every time, so it is a safe,
    always-present input edge).
    """
    input:
        feature_matrix=str(DOWNSTREAM_FEATURE_MATRIX),
        importances=str(DOWNSTREAM_MODELS_ROOT / "feature_importances_{feature_set}.json"),
    output:
        report=str(DOWNSTREAM_MODELS_ROOT / "tf_eval_report_{feature_set}.json"),
    params:
        model_root=lambda wc, output: str(Path(output.report).parent),
    log:
        str(DOWNSTREAM_MIXTURES_ROOT.parent / "models" / "tf_model_reshaped" / "logs" / "tf_eval_{feature_set}.log"),
    conda:
        "flash_fraq"
    threads: 2
    resources:
        mem_mb=8192,
        runtime=60,
    shell:
        f"mkdir -p $(dirname {{log}}) && {FLASHFRAQ_PYTHON} {SCRIPTS_DIR / 'run_tf_eval.py'} "
        "--feature-matrix-path {input.feature_matrix} "
        "--model-root {params.model_root} "
        "--report-path {output.report} "
        "--feature-set {wildcards.feature_set} "
        "> {log} 2>&1"


rule downstream_tf_eval_complete:
    input:
        expand(
            str(DOWNSTREAM_MODELS_ROOT / "tf_eval_report_{feature_set}.json"),
            feature_set=DOWNSTREAM_FEATURE_SETS,
        ),
    output:
        touch(str(DOWNSTREAM_MIXTURES_ROOT / "downstream_tf_eval_complete.flag")),
    conda:
        "flash_fraq"
