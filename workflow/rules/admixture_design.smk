"""Admixture design: native-coverage harvest, feasibility table, mixture plan.

Workstation-local, like cohort_qc: it reads the local_features derived tree and
writes four small CSVs. Nothing here needs HCC1395 or the cluster -- the tumor
source enters as config, not as an input file.

Ask for it by name:
  snakemake --workflow-profile workflow/profiles/local admixture_design_complete

Globals are prefixed ADMIXDES_ for the reason local_features.smk documents:
include: shares one Python namespace across every rule file.
"""

from pathlib import Path

ADMIXDES_CONFIG = config["admixture_design"]
ADMIXDES_DERIVED_ROOT = Path(ADMIXDES_CONFIG["derived_root"])
ADMIXDES_OUTPUT_ROOT = Path(ADMIXDES_CONFIG["output_root"])
ADMIXDES_COVERAGE_MANIFEST = ADMIXDES_OUTPUT_ROOT / "native_coverage.csv"
ADMIXDES_PLAN = ADMIXDES_OUTPUT_ROOT / "admixture_plan.csv"
ADMIXDES_FEASIBILITY = ADMIXDES_OUTPUT_ROOT / "admixture_feasibility.csv"
ADMIXDES_SPLITS = ADMIXDES_OUTPUT_ROOT / "split_assignment.csv"

ADMIXDES_DATASET_ARGUMENTS = " ".join(
    "--dataset {name}:{platform}:{source_group}:{tumor_status}".format(**dataset)
    for dataset in ADMIXDES_CONFIG["datasets"]
)
# Repo-relative in the config; resolved against REPO_ROOT because drivers run
# Snakemake from snakemake_workdirs/<name>, not from the repo.
ADMIXDES_PHENOTYPE_ARGUMENTS = " ".join(
    f"--phenotype {dataset['name']}={REPO_ROOT / dataset['phenotype_path']}"
    for dataset in ADMIXDES_CONFIG["datasets"]
    if dataset.get("phenotype_path")
)
ADMIXDES_PHENOTYPE_INPUTS = [
    str(REPO_ROOT / dataset["phenotype_path"]) for dataset in ADMIXDES_CONFIG["datasets"] if dataset.get("phenotype_path")
]
ADMIXDES_BACKGROUND_ARGUMENTS = " ".join(
    f"--background-source-group {group}" for group in ADMIXDES_CONFIG["background_source_groups"]
) + (
    f" --background-tumor-status {ADMIXDES_CONFIG['background_tumor_status']}"
    if ADMIXDES_CONFIG.get("background_tumor_status")
    else ""
)
ADMIXDES_SPLIT_ARGUMENTS = " ".join(
    f"--split-fraction {name}={weight}" for name, weight in ADMIXDES_CONFIG["split_fractions"].items()
)
ADMIXDES_LADDER = " ".join(str(value) for value in ADMIXDES_CONFIG["fraction_ladder"])
ADMIXDES_COVERAGES = " ".join(str(value) for value in ADMIXDES_CONFIG["candidate_target_coverages"])


rule native_coverage_manifest:
    """Genome-wide mean depth per sample, from the mosdepth summaries.

    Depends on local_features_complete rather than on the per-sample summaries
    for the same reason cohort_feature_matrix does: the harvester discovers
    samples by walking the derived tree, so naming its inputs individually would
    duplicate that discovery in the DAG and let the two disagree.
    """
    input:
        flag=str(LOCALFEAT_OUTPUT_ROOT / "local_features_complete.flag"),
        phenotypes=ADMIXDES_PHENOTYPE_INPUTS,
    output:
        manifest=str(ADMIXDES_COVERAGE_MANIFEST),
        summary=str(ADMIXDES_OUTPUT_ROOT / "native_coverage_summary.json"),
    log:
        str(ADMIXDES_OUTPUT_ROOT / "logs" / "native_coverage.log"),
    conda:
        "flash_fraq"
    threads: 1
    resources:
        mem_mb=2048,
        runtime=30,
    shell:
        f"python {SCRIPTS_DIR / 'run_native_coverage.py'} "
        f"--derived-root {ADMIXDES_DERIVED_ROOT} "
        f"--cram-root {ADMIXDES_CONFIG['cram_root']} "
        f"{ADMIXDES_DATASET_ARGUMENTS} "
        f"{ADMIXDES_PHENOTYPE_ARGUMENTS} "
        "--output-path {output.manifest} "
        "--summary-path {output.summary} "
        f"--reference-build {ADMIXDES_CONFIG['reference_build']} "
        "> {log} 2>&1"


checkpoint admixture_design:
    """Feasibility table, C_target choice, split assignment and mixture plan.

    Read the summary JSON, not just the CSVs: it carries the retention the
    chosen C_target buys at the lowest rung and the per-f source proportions,
    which are the two properties the design exists to guarantee.

    A checkpoint, not a plain rule: rule admixture_mixture (admixture.smk)
    needs one job per admixture_plan.csv row, but those rows -- and their
    count -- aren't known until this rule has actually run, so its consumer
    re-evaluates via checkpoints.admixture_design.get(...) rather than a
    statically expanded wildcard list.
    """
    input:
        manifest=str(ADMIXDES_COVERAGE_MANIFEST),
    output:
        feasibility=str(ADMIXDES_FEASIBILITY),
        plan=str(ADMIXDES_PLAN),
        splits=str(ADMIXDES_SPLITS),
        summary=str(ADMIXDES_OUTPUT_ROOT / "admixture_design_summary.json"),
    log:
        str(ADMIXDES_OUTPUT_ROOT / "logs" / "admixture_design.log"),
    conda:
        "flash_fraq"
    threads: 1
    resources:
        mem_mb=4096,
        runtime=30,
    shell:
        f"python {SCRIPTS_DIR / 'run_admixture_design.py'} "
        "--coverage-manifest-path {input.manifest} "
        f"--output-root {ADMIXDES_OUTPUT_ROOT} "
        "--summary-path {output.summary} "
        f"--tumor-sample-id {ADMIXDES_CONFIG['tumor']['sample_id']} "
        f"--tumor-subject-id {ADMIXDES_CONFIG['tumor']['subject_id']} "
        f"--tumor-platform {ADMIXDES_CONFIG['tumor']['platform']} "
        f"--tumor-native-coverage {ADMIXDES_CONFIG['tumor']['native_coverage']} "
        f"--tumor-bam-path {ADMIXDES_CONFIG['tumor']['bam_path']} "
        f"--matched-pair-sample-id {ADMIXDES_CONFIG['matched_pair_background']['sample_id']} "
        f"--matched-pair-subject-id {ADMIXDES_CONFIG['matched_pair_background']['subject_id']} "
        f"--matched-pair-source-group {ADMIXDES_CONFIG['matched_pair_background']['source_group']} "
        f"--matched-pair-platform {ADMIXDES_CONFIG['matched_pair_background']['platform']} "
        f"--matched-pair-native-coverage {ADMIXDES_CONFIG['matched_pair_background']['native_coverage']} "
        f"--matched-pair-bam-path {ADMIXDES_CONFIG['matched_pair_background']['bam_path']} "
        f"--matched-pair-replicates-per-fraction {ADMIXDES_CONFIG['matched_pair_replicates_per_fraction']} "
        f"--fraction-ladder {ADMIXDES_LADDER} "
        f"--candidate-target-coverages {ADMIXDES_COVERAGES} "
        f"{ADMIXDES_BACKGROUND_ARGUMENTS} "
        f"{ADMIXDES_SPLIT_ARGUMENTS} "
        f"--mixtures-per-fraction {ADMIXDES_CONFIG['mixtures_per_fraction']} "
        f"--min-background-samples-per-source {ADMIXDES_CONFIG['min_background_samples_per_source']} "
        f"--min-background-retention-per-source {ADMIXDES_CONFIG['min_background_retention_per_source']} "
        f"--reference-build {ADMIXDES_CONFIG['reference_build']} "
        f"--seed {ADMIXDES_CONFIG['seed']} "
        "> {log} 2>&1"


rule admixture_design_complete:
    input:
        manifest=str(ADMIXDES_COVERAGE_MANIFEST),
        feasibility=str(ADMIXDES_FEASIBILITY),
        plan=str(ADMIXDES_PLAN),
        splits=str(ADMIXDES_SPLITS),
    output:
        touch(str(ADMIXDES_OUTPUT_ROOT / "admixture_design_complete.flag")),
    conda:
        "flash_fraq"
