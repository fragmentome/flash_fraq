#!/usr/bin/env bash
# Drive the full downstream rerun on the reshaped HCC1395 admixture mixtures
# (docs/todo/todo/reshaped-mixture-downstream-rerun.md), once
# reshaped_tumor_and_mixture_chain.sh's live mixture build finishes.
#
# Waits on that chain's own log for its exact completion line, verifies every
# real_background mixture has the inputs downstream needs, then runs, in
# ticket order and each independently resumable:
#   relabel -> fragments -> gc_correct -> ichorCNA(low_tf) -> feature_matrix
#   -> decoy_qc (HARD gate) -> tf_model(all, cnv_only, fragment_only)
#   -> tf_eval(all, cnv_only, fragment_only)
#
# The four per-mixture stages reuse the existing, already-tested bash+xargs
# drivers (admixture_mixture_{relabel,fragments,gc_correct,ichor}_driver.sh),
# pointed at the reshaped mixture root via their MIXTURES_ROOT/STAGING_ROOT env
# overrides -- see each driver's own docstring for why they are not Snakemake
# rules (flat mix_<name>/ layout, not the cohort CRAM-mirror layout
# workflow/rules/local_features.smk expects). The four whole-matrix stages run
# through a brand-new, separate Snakemake entry point,
# workflow/Snakefile_reshaped_downstream (+ workflow/rules/admixture_downstream.smk
# + workflow/config/admixture_downstream_reshaped.yaml), driven from its own
# workdir (snakemake_workdirs/reshaped_downstream) -- this NEVER shares a DAG,
# workdir, or config with the live snakemake_workdirs/admixture build driven by
# the main workflow/Snakefile; this script never invokes snakemake against that
# workdir, and never kills or otherwise touches that build's process.
#
# Everything here only ever READS the reshaped mixture build's outputs and the
# original (unreshaped) cohort's matched_pair plan rows, and only ever WRITES
# under admixture_mixtures_reshaped/ and models/tf_model_reshaped/ -- distinct
# output roots from the original 495-mixture run. Nothing is ever deleted by
# this script.
#
# Hard gate: decoy QC must pass (decoy_qc_report_real_background.json's
# `passed` field, per src/flash_fraq/decoy_qc.py's own report_path()/`passed`
# semantics) before tf_model/tf_eval run. If it does not pass, this script
# stops with a clear message and a non-zero exit, same as the original
# 495-mixture run's own decoy-QC gate.
#
# Env hooks:
#   SKIP_WAIT=1        skip the wait-gate on the mixture-build log entirely
#                      (e.g. to test this script, or if the build already
#                      finished and its log is gone/rotated).
#   DRY_RUN=1          the four Snakemake-backed stages (feature_matrix,
#                      decoy_qc, tf_model, tf_eval) run with `snakemake -n`
#                      only, printing each stage's job stats; nothing is
#                      executed. The four per-mixture bash-driver stages are
#                      NOT invoked at all under DRY_RUN=1 (they are plain
#                      bash+xargs loops, not Snakemake rules, so there is no
#                      `-n` equivalent for them) -- this script instead prints
#                      how many mixture directories exist and how many already
#                      look done for each of those stages, and skips the
#                      real_background-completeness gate's hard-abort (it
#                      still reports what it finds, since not every mixture
#                      will exist yet while the upstream build is mid-flight).
#   STOP_AFTER=<step>  stop (exit 0) right after the named step completes.
#                      One of: wait, verify, relabel, fragments, gc_correct,
#                      ichor, combined_plan, feature_matrix, decoy_qc,
#                      tf_model, tf_eval.
#   OVERRIDE_DECOY_GATE=1  log the decoy-QC verdict but continue to
#                      tf_model/tf_eval even when `passed` is false. For a
#                      gate shown to be miscalibrated (2026-10-04: the
#                      library-confound check exceeds |0.1| on all 50
#                      label permutations, null rho -0.27 +/- 0.04, observed
#                      -0.117); results from such a run must be reported as
#                      gate-overridden.
#
# Usage (wait for the live mixture build, then run everything):
#   setsid nohup ./scripts/drivers/reshaped_downstream_chain.sh \
#       > logs/reshaped_downstream_chain_nohup_$(date +%Y%m%d_%H%M%S).log 2>&1 &
#
# Testing (never against the live build/workdir):
#   bash -n scripts/drivers/reshaped_downstream_chain.sh
#   DRY_RUN=1 SKIP_WAIT=1 ./scripts/drivers/reshaped_downstream_chain.sh

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"

ts() { date -Is; }

mkdir -p logs
LOG="logs/reshaped_downstream_chain_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1

trap 'echo "[$(ts)] reshaped_downstream_chain CHAIN FAILED at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR

echo "[$(ts)] reshaped_downstream_chain started; log ${LOG}"

# ---------------------------------------------------------------------------
# Configuration (all env-overridable; defaults match the live reshaped build).
# ---------------------------------------------------------------------------

# Wait-gate on the upstream mixture-build chain.
MIXTURE_CHAIN_LOG="${MIXTURE_CHAIN_LOG:-logs/reshaped_tumor_and_mixture_chain_nohup_20261001_191542.log}"
WAIT_POLL_SECONDS="${WAIT_POLL_SECONDS:-60}"
SKIP_WAIT="${SKIP_WAIT:-0}"

DRY_RUN="${DRY_RUN:-0}"
STOP_AFTER="${STOP_AFTER:-}"
OVERRIDE_DECOY_GATE="${OVERRIDE_DECOY_GATE:-0}"

# Reshaped mixture build outputs (read) / downstream outputs (written).
MIXTURES_ROOT="${MIXTURES_ROOT:-/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures_reshaped}"
COHORT_RESHAPED_ROOT="${COHORT_RESHAPED_ROOT:-/mnt/8tb_hdd/flash_fraq_derived/cohort_reshaped}"
PLAN_PATH="${PLAN_PATH:-${COHORT_RESHAPED_ROOT}/admixture_plan.csv}"
SPLIT_ASSIGNMENT_PATH="${SPLIT_ASSIGNMENT_PATH:-${COHORT_RESHAPED_ROOT}/split_assignment.csv}"
COMBINED_PLAN_PATH="${COMBINED_PLAN_PATH:-${COHORT_RESHAPED_ROOT}/admixture_plan_combined.csv}"

# Original (unreshaped) cohort's plan -- read-only, source of the
# matched_pair-only rows the combined plan needs (see below).
ORIGINAL_PLAN_PATH="${ORIGINAL_PLAN_PATH:-/mnt/8tb_hdd/flash_fraq_derived/cohort/admixture_plan.csv}"

# Per-mixture bash-driver staging roots (unused directly here; passed through
# to the drivers themselves, which already default sensibly).

# Snakemake-backed whole-matrix stages.
DOWNSTREAM_SNAKEFILE="${REPO_ROOT}/workflow/Snakefile_reshaped_downstream"
DOWNSTREAM_WORKDIR="${REPO_ROOT}/snakemake_workdirs/reshaped_downstream"
SNAKEMAKE_CORES="${SNAKEMAKE_CORES:-36}"
SNAKEMAKE_MEM_MB="${SNAKEMAKE_MEM_MB:-40000}"

# Per-mixture bash-driver concurrency (each driver also has its own default;
# these are passed through only so one place documents/overrides all of them).
RELABEL_JOBS="${RELABEL_JOBS:-8}"
FRAGMENTS_JOBS="${FRAGMENTS_JOBS:-8}"
GC_CORRECT_JOBS="${GC_CORRECT_JOBS:-16}"
ICHOR_JOBS="${ICHOR_JOBS:-6}"
ICHOR_PROFILE="${ICHOR_PROFILE:-low_tf}"

step_done() {
    # Print a uniform "finished step X" line and honor STOP_AFTER.
    local step="$1"
    echo "[$(ts)] reshaped_downstream_chain: step '${step}' complete"
    if [ "${STOP_AFTER}" = "${step}" ]; then
        echo "[$(ts)] reshaped_downstream_chain: STOP_AFTER=${step} -- stopping here"
        exit 0
    fi
}

# ---------------------------------------------------------------------------
# Step: wait
#
# Poll MIXTURE_CHAIN_LOG for the upstream chain's own exact completion line.
# Must match "reshaped_tumor_and_mixture_chain CHAIN DONE" exactly -- the same
# log also contains an unrelated line from that chain's own wait-gate on an
# *earlier* run, "CHAIN DONE observed in .../upstream_chain_*.log; proceeding",
# which must NOT match (it does not contain the substring
# "reshaped_tumor_and_mixture_chain CHAIN DONE" at all, so a plain fixed-string
# grep is safe and does not need anchoring).
# ---------------------------------------------------------------------------
if [ "${SKIP_WAIT}" = "1" ]; then
    echo "[$(ts)] reshaped_downstream_chain: SKIP_WAIT=1 -- not waiting on ${MIXTURE_CHAIN_LOG}"
else
    echo "[$(ts)] reshaped_downstream_chain: waiting for 'reshaped_tumor_and_mixture_chain CHAIN DONE' in ${MIXTURE_CHAIN_LOG} (polling every ${WAIT_POLL_SECONDS}s)"
    while true; do
        if [ -s "${MIXTURE_CHAIN_LOG}" ]; then
            if grep -qF 'reshaped_tumor_and_mixture_chain CHAIN FAILED' "${MIXTURE_CHAIN_LOG}"; then
                echo "[$(ts)] ABORT: ${MIXTURE_CHAIN_LOG} reports 'reshaped_tumor_and_mixture_chain CHAIN FAILED' -- not proceeding" >&2
                exit 1
            fi
            if grep -qF 'reshaped_tumor_and_mixture_chain CHAIN DONE' "${MIXTURE_CHAIN_LOG}"; then
                echo "[$(ts)] reshaped_downstream_chain: observed CHAIN DONE in ${MIXTURE_CHAIN_LOG}; proceeding"
                break
            fi
        fi
        sleep "${WAIT_POLL_SECONDS}"
    done
fi
step_done wait

# ---------------------------------------------------------------------------
# Step: verify
#
# Every real_background mix_* dir named in PLAN_PATH must have its dedup BAM
# (+ index) and sample_metadata.csv before any downstream stage runs. Under
# DRY_RUN=1, missing mixtures are reported but do not abort the script (the
# mixtures may genuinely not all exist yet while the upstream build is
# mid-flight; a dry run against the existing subset is still useful).
# ---------------------------------------------------------------------------
echo "[$(ts)] reshaped_downstream_chain: verifying real_background mixtures under ${MIXTURES_ROOT} against ${PLAN_PATH}"

[ -s "${PLAN_PATH}" ] || { echo "[$(ts)] ABORT: missing ${PLAN_PATH}" >&2; exit 1; }
[ -d "${MIXTURES_ROOT}" ] || { echo "[$(ts)] ABORT: missing ${MIXTURES_ROOT}" >&2; exit 1; }

VERIFY_STATUS=0
"${FLASHFRAQ_PYTHON:-/home/matvii/miniconda3/envs/flash_fraq/bin/python}" - "${PLAN_PATH}" "${MIXTURES_ROOT}" <<'PYEOF' || VERIFY_STATUS=$?
import csv
import sys
from pathlib import Path

plan_path, mixtures_root = sys.argv[1], Path(sys.argv[2])
with open(plan_path) as f:
    rows = [r for r in csv.DictReader(f) if r.get("regime") == "real_background"]

missing = []
for row in rows:
    sample = row["sample_name"]
    mix_dir = mixtures_root / sample
    dedup_bam = mix_dir / f"{sample}.dedup.bam"
    metadata = mix_dir / "sample_metadata.csv"
    problems = []
    if not dedup_bam.is_file() or dedup_bam.stat().st_size == 0:
        problems.append("missing/empty dedup BAM")
    if not (dedup_bam.with_suffix(dedup_bam.suffix + ".bai")).is_file():
        problems.append("missing .bai index")
    if not metadata.is_file() or metadata.stat().st_size == 0:
        problems.append("missing/empty sample_metadata.csv")
    if problems:
        missing.append((sample, problems))

print(f"  {len(rows)} real_background plan rows, {len(missing)} with missing outputs")
if missing:
    for sample, problems in missing[:50]:
        print(f"  MISSING {sample}: {', '.join(problems)}")
    if len(missing) > 50:
        print(f"  ... and {len(missing) - 50} more")
    sys.exit(1)
PYEOF

if [ "${VERIFY_STATUS}" -ne 0 ]; then
    if [ "${DRY_RUN}" = "1" ]; then
        echo "[$(ts)] reshaped_downstream_chain: DRY_RUN=1 -- some real_background mixtures are incomplete (expected mid-build); continuing with dry runs against whatever exists"
    else
        echo "[$(ts)] ABORT: real_background mixtures are incomplete -- see MISSING lines above. Re-run once reshaped_tumor_and_mixture_chain.sh has truly finished." >&2
        exit 1
    fi
fi
step_done verify

# ---------------------------------------------------------------------------
# Step: combined_plan
#
# scripts/run_admixture_mixture_feature_matrix.py aborts if any mix_* sample
# has no matching row in --admixture-plan-path. The reshaped PLAN_PATH has
# zero matched_pair rows (matched_pair_replicates_per_fraction: 0 in
# admixture_design_reshaped.yaml -- those 11 mix_* dirs under MIXTURES_ROOT are
# symlinks to the ORIGINAL, unreshaped build), so a plan built from PLAN_PATH
# alone would fail on them. Build a combined plan = reshaped plan's
# real_background rows + the original cohort's matched_pair-only rows. This is
# cheap, so it is rebuilt on every invocation (idempotent, never a
# skip-if-done gate).
# ---------------------------------------------------------------------------
echo "[$(ts)] reshaped_downstream_chain: building combined plan at ${COMBINED_PLAN_PATH}"
[ -s "${ORIGINAL_PLAN_PATH}" ] || { echo "[$(ts)] ABORT: missing ${ORIGINAL_PLAN_PATH}" >&2; exit 1; }

"${FLASHFRAQ_PYTHON:-/home/matvii/miniconda3/envs/flash_fraq/bin/python}" - "${PLAN_PATH}" "${ORIGINAL_PLAN_PATH}" "${COMBINED_PLAN_PATH}" <<'PYEOF'
import csv
import sys

reshaped_plan_path, original_plan_path, combined_plan_path = sys.argv[1:4]

with open(reshaped_plan_path) as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    reshaped_rows = [r for r in reader if r.get("regime") == "real_background"]

with open(original_plan_path) as f:
    original_rows = [r for r in csv.DictReader(f) if r.get("regime") == "matched_pair"]

combined_rows = reshaped_rows + original_rows
with open(combined_plan_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for row in combined_rows:
        writer.writerow({k: row.get(k, "") for k in fieldnames})

print(
    f"  wrote {len(combined_rows)} rows to {combined_plan_path} "
    f"({len(reshaped_rows)} real_background from {reshaped_plan_path}, "
    f"{len(original_rows)} matched_pair from {original_plan_path})"
)
PYEOF
step_done combined_plan

# ---------------------------------------------------------------------------
# Per-mixture bash-driver stages: relabel -> fragments -> gc_correct -> ichor.
#
# Each driver is already idempotent (skip-if-done) and SSD-stages its own BAM
# I/O; none of them are re-implemented here. Under DRY_RUN=1 they are not
# invoked (see header comment) -- instead this prints how many mix_* dirs
# exist and how many already look done for each stage, so DRY_RUN still gives
# a useful signal without risking hours of real work.
# ---------------------------------------------------------------------------

dry_run_stage_report() {
    local label="$1" done_glob="$2"
    local total done
    total=$(find "${MIXTURES_ROOT}" -mindepth 1 -maxdepth 1 -type d -name 'mix_*' | wc -l)
    done=$(find "${MIXTURES_ROOT}" -mindepth 2 -maxdepth 3 -name "${done_glob}" 2>/dev/null | wc -l)
    echo "[$(ts)] reshaped_downstream_chain: DRY_RUN ${label}: ${done}/${total} mix_* dirs already look done (glob: ${done_glob})"
}

if [ "${DRY_RUN}" = "1" ]; then
    echo "[$(ts)] reshaped_downstream_chain: DRY_RUN=1 -- skipping real relabel/fragments/gc_correct/ichor execution (plain bash+xargs drivers, not Snakemake rules; no -n equivalent)"
    dry_run_stage_report relabel '*.empirical_tumor_fraction.sites.sites.tsv.gz'
    step_done relabel
    dry_run_stage_report fragments 'fragment_stats.csv'
    step_done fragments
    dry_run_stage_report gc_correct 'arm_zscores.csv'
    step_done gc_correct
    dry_run_stage_report ichor 'ichor_estimates.csv'
    step_done ichor
else
    echo "[$(ts)] reshaped_downstream_chain: relabel"
    MIXTURES_ROOT="${MIXTURES_ROOT}" JOBS="${RELABEL_JOBS}" \
        bash "${REPO_ROOT}/scripts/drivers/admixture_mixture_relabel_driver.sh"
    step_done relabel

    echo "[$(ts)] reshaped_downstream_chain: fragments"
    MIXTURES_ROOT="${MIXTURES_ROOT}" JOBS="${FRAGMENTS_JOBS}" \
        bash "${REPO_ROOT}/scripts/drivers/admixture_mixture_fragments_driver.sh"
    step_done fragments

    echo "[$(ts)] reshaped_downstream_chain: gc_correct"
    MIXTURES_ROOT="${MIXTURES_ROOT}" JOBS="${GC_CORRECT_JOBS}" \
        bash "${REPO_ROOT}/scripts/drivers/admixture_mixture_gc_correct_driver.sh"
    step_done gc_correct

    echo "[$(ts)] reshaped_downstream_chain: ichor (profile=${ICHOR_PROFILE})"
    MIXTURES_ROOT="${MIXTURES_ROOT}" JOBS="${ICHOR_JOBS}" ICHOR_PROFILE="${ICHOR_PROFILE}" \
        bash "${REPO_ROOT}/scripts/drivers/admixture_mixture_ichor_driver.sh"
    step_done ichor
fi

# ---------------------------------------------------------------------------
# Snakemake-backed whole-matrix stages: feature_matrix -> decoy_qc (HARD gate)
# -> tf_model -> tf_eval, via the standalone workflow/Snakefile_reshaped_downstream.
#
# Deliberately a separate workdir/Snakefile from the live
# snakemake_workdirs/admixture build -- never touched here.
# ---------------------------------------------------------------------------

. "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate snakemake
export PATH="${HOME}/miniconda3/envs/mosdepth/bin:${PATH}"
export PYTHONPATH="${REPO_ROOT}/src"

mkdir -p "${DOWNSTREAM_WORKDIR}"

echo "[$(ts)] reshaped_downstream_chain: snakemake --lint on ${DOWNSTREAM_SNAKEFILE}"
snakemake --lint -s "${DOWNSTREAM_SNAKEFILE}"

SNAKEMAKE_ARGS=(-p -s "${DOWNSTREAM_SNAKEFILE}"
                --cores "${SNAKEMAKE_CORES}"
                --resources "mem_mb=${SNAKEMAKE_MEM_MB}"
                --rerun-triggers mtime params input code)

run_snakemake_target() {
    # $1 = step name (for STOP_AFTER/logging), $2.. = snakemake target(s).
    local step="$1"; shift
    local targets=("$@")
    (
        cd "${DOWNSTREAM_WORKDIR}"
        if [ "${DRY_RUN}" = "1" ]; then
            echo "[$(ts)] reshaped_downstream_chain: DRY_RUN dry-run for step '${step}' (targets: ${targets[*]})"
            snakemake -n "${SNAKEMAKE_ARGS[@]}" -- "${targets[@]}"
        else
            snakemake "${SNAKEMAKE_ARGS[@]}" -- "${targets[@]}"
        fi
    )
}

echo "[$(ts)] reshaped_downstream_chain: feature_matrix"
# The feature-matrix rule's output is declared by path in
# workflow/config/admixture_downstream_reshaped.yaml; target it directly.
run_snakemake_target feature_matrix \
    "/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures_reshaped/feature_matrix.csv"
step_done feature_matrix

echo "[$(ts)] reshaped_downstream_chain: decoy_qc (both regimes)"
run_snakemake_target decoy_qc \
    "/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures_reshaped/decoy_qc/decoy_qc_report_real_background.json" \
    "/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures_reshaped/decoy_qc/decoy_qc_report_matched_pair.json"
step_done decoy_qc

if [ "${DRY_RUN}" = "1" ]; then
    echo "[$(ts)] reshaped_downstream_chain: DRY_RUN=1 -- skipping the decoy-QC hard-gate check (report may not exist yet) and continuing to dry-run tf_model/tf_eval"
else
    DECOY_QC_REPORT="/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures_reshaped/decoy_qc/decoy_qc_report_real_background.json"
    echo "[$(ts)] reshaped_downstream_chain: checking HARD gate -- ${DECOY_QC_REPORT}'s passed field"
    GATE_PASSED="$("${FLASHFRAQ_PYTHON:-/home/matvii/miniconda3/envs/flash_fraq/bin/python}" -c "
import json
with open('${DECOY_QC_REPORT}') as f:
    report = json.load(f)
print('true' if report.get('passed') else 'false')
")"
    if [ "${GATE_PASSED}" != "true" ] && [ "${OVERRIDE_DECOY_GATE}" = "1" ]; then
        echo "[$(ts)] reshaped_downstream_chain: WARNING decoy QC did not pass for real_background (see ${DECOY_QC_REPORT}) -- OVERRIDE_DECOY_GATE=1, continuing; report results as gate-overridden"
    elif [ "${GATE_PASSED}" != "true" ]; then
        echo "[$(ts)] ABORT: decoy QC did not pass for real_background (see ${DECOY_QC_REPORT}) -- stopping before tf_model/tf_eval, same hard-gate semantics as src/flash_fraq/decoy_qc.py's run()/passed field." >&2
        exit 1
    else
        echo "[$(ts)] reshaped_downstream_chain: decoy QC passed for real_background -- proceeding"
    fi
fi

echo "[$(ts)] reshaped_downstream_chain: tf_model (feature sets: all, cnv_only, fragment_only)"
run_snakemake_target tf_model \
    "/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_reshaped/feature_importances_all.json" \
    "/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_reshaped/feature_importances_cnv_only.json" \
    "/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_reshaped/feature_importances_fragment_only.json"
step_done tf_model

echo "[$(ts)] reshaped_downstream_chain: tf_eval (feature sets: all, cnv_only, fragment_only)"
run_snakemake_target tf_eval \
    "/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_reshaped/tf_eval_report_all.json" \
    "/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_reshaped/tf_eval_report_cnv_only.json" \
    "/mnt/8tb_hdd/flash_fraq_derived/models/tf_model_reshaped/tf_eval_report_fragment_only.json"
step_done tf_eval

echo "[$(ts)] reshaped_downstream_chain CHAIN DONE"
