"""Cohort feature matrix + phase2 §6 decoy-classifier batch-effect gate.

Runs after local_features, on the same workstation and over the same derived
tree: it consumes nothing but the per-sample tables that workflow already wrote,
so there is no staging, no CRAM access and no cluster dependency here.

Ask for it by name, like local_features (the Snakefile's `rule all` covers the
cluster-side alignment targets only):
  snakemake --workflow-profile workflow/profiles/local cohort_qc_complete

Globals are prefixed COHORTQC_ for the reason local_features.smk documents:
include: shares one Python namespace across every rule file.
"""

from pathlib import Path

COHORTQC_CONFIG = config["cohort_qc"]
COHORTQC_DERIVED_ROOT = Path(COHORTQC_CONFIG["derived_root"])
COHORTQC_OUTPUT_ROOT = Path(COHORTQC_CONFIG["output_root"])
COHORTQC_DATASET_PLATFORMS = COHORTQC_CONFIG["dataset_platforms"]
COHORTQC_MATRIX = COHORTQC_OUTPUT_ROOT / "cohort_feature_matrix.csv"
COHORTQC_SUMMARY = COHORTQC_OUTPUT_ROOT / "cohort_feature_matrix_summary.json"
COHORTQC_REPORT = COHORTQC_OUTPUT_ROOT / "decoy_qc" / "decoy_qc_report.json"

COHORTQC_DATASET_ARGUMENTS = " ".join(
    f"--dataset {dataset}={platform}" for dataset, platform in COHORTQC_DATASET_PLATFORMS.items()
)


rule cohort_feature_matrix:
    """One row per real (unmixed) sample, joined across MAT-14/15/16 outputs.

    Depends on local_features_complete rather than on the individual per-sample
    tables: the assembler discovers samples by walking the derived tree, so
    naming its inputs one by one would duplicate that discovery in the DAG and
    let the two disagree. The cohort-level flag is the honest dependency -- the
    matrix is only meaningful once every sample's tables are written.
    """
    input:
        flag=str(LOCALFEAT_OUTPUT_ROOT / "local_features_complete.flag"),
    output:
        matrix=str(COHORTQC_MATRIX),
        summary=str(COHORTQC_SUMMARY),
    log:
        str(COHORTQC_OUTPUT_ROOT / "logs" / "cohort_feature_matrix.log"),
    conda:
        "flash_fraq"
    threads: 1
    resources:
        mem_mb=16384,
        runtime=60,
    shell:
        f"python {SCRIPTS_DIR / 'run_cohort_feature_matrix.py'} "
        f"--derived-root {COHORTQC_DERIVED_ROOT} "
        f"{COHORTQC_DATASET_ARGUMENTS} "
        "--output-path {output.matrix} "
        "--summary-path {output.summary} "
        f"--reference-build {COHORTQC_CONFIG['reference_build']} "
        "> {log} 2>&1"


rule cohort_decoy_qc:
    """phase2 §6: can source/platform be predicted from the features, and does
    that prediction correlate with the tumor-fraction label?

    §6 says to run this "before any real model training is treated as
    informative". At this stage the label is ichorCNA's estimate rather than an
    SNV-VAF empirical fraction, which makes it a weaker check than the one Phase
    3 needs -- hence fail_on_flag defaulting to false in the config.
    """
    input:
        matrix=str(COHORTQC_MATRIX),
    output:
        report=str(COHORTQC_REPORT),
    params:
        output_dir=lambda wc, output: str(Path(output.report).parent),
        fail_on_flag="" if COHORTQC_CONFIG["fail_on_flag"] else "--no-fail-on-flag ",
    log:
        str(COHORTQC_OUTPUT_ROOT / "logs" / "cohort_decoy_qc.log"),
    conda:
        "flash_fraq"
    threads: 1
    resources:
        mem_mb=16384,
        runtime=60,
    shell:
        f"python {SCRIPTS_DIR / 'run_decoy_qc.py'} "
        "--feature-matrix-path {input.matrix} "
        "--output-dir {params.output_dir} "
        "--matrix-kind cohort "
        f"--correlation-threshold {COHORTQC_CONFIG['correlation_threshold']} "
        f"--fixed-fraction-target {COHORTQC_CONFIG['fixed_fraction_target']} "
        f"--fixed-fraction-tolerance {COHORTQC_CONFIG['fixed_fraction_tolerance']} "
        f"--max-missing-fraction {COHORTQC_CONFIG['max_missing_fraction']} "
        f"--n-splits {COHORTQC_CONFIG['n_splits']} "
        "{params.fail_on_flag}"
        "> {log} 2>&1"


rule cohort_qc_complete:
    input:
        matrix=str(COHORTQC_MATRIX),
        summary=str(COHORTQC_SUMMARY),
        report=str(COHORTQC_REPORT),
    output:
        touch(str(COHORTQC_OUTPUT_ROOT / "cohort_qc_complete.flag")),
    conda:
        "flash_fraq"
