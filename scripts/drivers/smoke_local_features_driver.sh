#!/usr/bin/env bash
# Smoke-test driver for MAT-14/15/16: run ONE sample through the full local_features
# chain (docs/local-features-pipeline.md) -- stage_cram -> verify_cram ->
# coverage/fragments/gc_correct -> cram_to_bam -> ichor -> sample_complete --
# on the workstation, before pointing local_features_complete at the whole
# Snyder+Cristiano cohort.
#
# Unlike scripts/drivers/mat19_driver.sh etc., this never SSHes anywhere --
# local_features is workstation-only (docs/local-features-pipeline.md §1: "the
# cluster's job ends at FASTQ->BAM alignment"), so this runs snakemake
# in-place against the local checkout.
#
# Default sample: Snyder_Shendure_WPS/IA04 -- smallest-looking CRAM (~1.2 GB)
# next to IA01, the sample docs/local-features-pipeline.md already validated
# end-to-end, so this exercises a second sample through the real Snakemake DAG
# (IA01's validation run went through the standalone scripts/run_local_*.py
# drivers directly, not `snakemake ... local_features_complete`).
#
# Usage:
#   ./scripts/drivers/smoke_local_features_driver.sh
#
# Overridable via env var: SMOKE_DATASET, SMOKE_SAMPLE.
#
# LOCALFEAT_OUTPUT_ROOT is hardcoded to match workflow/config/local_features.yaml's
# output_root, since the target path is built from it.
#
# Runs the CLAUDE.md-mandated pre-flight checks (`snakemake --lint`, then a
# dry run) before the real run. Requires the `snakemake` conda env (snakemake
# itself is not in the `flash_fraq` env -- see docs/local-features-pipeline.md
# §8) and PYTHONPATH=$PWD/src so the rule scripts can import flash_fraq.

set -euo pipefail

SMOKE_DATASET="${SMOKE_DATASET:-Snyder_Shendure_WPS}"
SMOKE_SAMPLE="${SMOKE_SAMPLE:-IA04}"

LOCALFEAT_OUTPUT_ROOT="/mnt/8tb_hdd/flash_fraq_derived"
TARGET="${LOCALFEAT_OUTPUT_ROOT}/${SMOKE_DATASET}/${SMOKE_SAMPLE}/localfeat_sample_complete.flag"

cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"
SNAKEFILE="${REPO_ROOT}/workflow/Snakefile"
LOCAL_PROFILE="${REPO_ROOT}/workflow/profiles/local"

. "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate snakemake
export PYTHONPATH="${REPO_ROOT}/src"

# Run from this workflow's own directory so its .snakemake -- lock, metadata,
# executor logs -- is separate from the cluster workflows' (see
# snakemake_workdirs/README.md). Paths above are absolute for that reason.
WORKFLOW_WORKDIR="${REPO_ROOT}/snakemake_workdirs/local_features"
mkdir -p "${WORKFLOW_WORKDIR}"
cd "${WORKFLOW_WORKDIR}"

echo "Pre-flight checks..."
snakemake --lint -s "${SNAKEFILE}"
snakemake -n -p --use-conda -s "${SNAKEFILE}" \
  --workflow-profile "${LOCAL_PROFILE}" "${TARGET}"

echo
echo "Running local_features for ${SMOKE_DATASET}/${SMOKE_SAMPLE}..."
time snakemake -p --use-conda -s "${SNAKEFILE}" \
  --workflow-profile "${LOCAL_PROFILE}" "${TARGET}"

echo
echo "Done. Expected artifacts under ${LOCALFEAT_OUTPUT_ROOT}/${SMOKE_DATASET}/${SMOKE_SAMPLE}/:"
echo "  cram_header_verification.json"
echo "  coverage/${SMOKE_SAMPLE}.mosdepth.summary.txt"
echo "  fragment_stats.csv, window_ratios.csv, arm_zscores.csv"
echo "  ichor/cnv_bins.csv, ichor/ichor_estimates.csv"
echo "  localfeat_sample_complete.flag"
