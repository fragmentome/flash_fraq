#!/usr/bin/env bash
# Drive the real in-silico admixture build: all 495 admixture_plan.csv rows
# (11 matched_pair + 484 real_background, MAT-13) through rule
# admixture_mixture (workflow/rules/admixture.smk) against the 18 real
# merged HCC1395/HCC1395BL library CRAMs and the Cristiano DELFI background
# pool -- this is the "real run against merged CRAMs" MAT-13's ticket has
# been waiting on since the 2026-09-17 checkpoint.
#
# Workstation-local, single Snakemake invocation, same shape as
# scripts/drivers/smoke_hcc1395_local_driver.sh: own snakemake_workdirs/
# subdirectory, `snakemake` conda env (not `flash_fraq`) for the CLI,
# PYTHONPATH=src so the rule files' `import flash_fraq...` resolve.
#
# NOT run through --workflow-profile workflow/profiles/local: this
# snakemake (9.11.4) throws `ValueError: not enough values to unpack` in
# resources.py's parse_resources when a CLI --resources list is combined
# with a profile config.yaml that also declares a `resources:` mapping (the
# local profile's ssd_mb/ena_download/sra_fetch/bam_align caps, none of
# which admixture_mixture consumes anyway). --cores/--resources/--keep-going
# are passed directly instead, which is a clean rule-level equivalent for a
# workflow that only touches one resource family (mem_mb) the profile
# doesn't model.
#
# Concurrency: each admixture_mixture job wants threads=8, mem_mb=8192
# (workflow/config/admixture.yaml). SNAKEMAKE_CORES=32 on a 36-core/125GB
# box caps concurrency at 4 jobs by the thread budget alone (32/8=4).
# SNAKEMAKE_MEM_MB must be >= 4*8192=32768 or the mem_mb resource cap binds
# first and silently drops concurrency to 3 -- caught live on the first
# launch (32000 fit only 3 jobs), hence the default below has headroom.
#
# Resilience: like hcc1395_vacation_queue_launcher.sh, this runs several
# passes rather than one. --keep-going means one bad mixture row does not
# abort the other 494; a later pass only re-attempts rows whose outputs are
# still missing (Snakemake's own file-based tracking is idempotent here --
# there is no separate already_done() to write, unlike the ENA fetch driver).
#
# Usage:
#   nohup ./scripts/drivers/admixture_mixtures_local_driver.sh \
#       > logs/admixture_mixtures_YYYYMMDD_HHMMSS.log 2>&1 &
#   (the script also tees its own timestamped log under logs/, so redirecting
#   stdout is optional -- see LOG below.)
#
# Env: SNAKEMAKE_CORES (32), SNAKEMAKE_MEM_MB (40000), PASSES (3),
# MIN_FREE_GB (100 -- generous given output_root's ~5-10G peak scratch per
# mixture against 3+ TB free on /mnt/8tb_hdd).
#
# Reused as-is for the reshaped-tumor rebuild
# (docs/todo/todo/reshaped-mixture-downstream-rerun.md,
# scripts/drivers/reshaped_tumor_and_mixture_chain.sh): ADMIX_OUTPUT_ROOT,
# PLAN_PATH, QUICKCHECK_CRAMS and EXTRA_CONFIGFILES are now env-overridable
# (defaults below reproduce the original, unparameterized behaviour exactly)
# rather than forking a second copy of this driver for a different
# output_root/plan/tumor source. QUICKCHECK_CRAMS is the pre-flight sanity
# list -- the reshaped run points it at the reshaped tumor BAM instead of the
# raw NV_3 CRAM, since that is this run's actual tumor source; the *generic*
# per-row path check a few lines down is dataset-agnostic and needs no
# override either way.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"

SNAKEFILE="${REPO_ROOT}/workflow/Snakefile"
ADMIX_OUTPUT_ROOT="${ADMIX_OUTPUT_ROOT:-/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures}"
PLAN_PATH="${PLAN_PATH:-/mnt/8tb_hdd/flash_fraq_derived/cohort/admixture_plan.csv}"
TARGET="${ADMIX_OUTPUT_ROOT}/admixture_mixtures_complete.flag"

SNAKEMAKE_CORES="${SNAKEMAKE_CORES:-32}"
SNAKEMAKE_MEM_MB="${SNAKEMAKE_MEM_MB:-40000}"
PASSES="${PASSES:-3}"
MIN_FREE_GB="${MIN_FREE_GB:-100}"
# Space-separated list of CRAM/BAM paths for the pre-flight quickcheck below.
# shellcheck disable=SC2206
QUICKCHECK_CRAMS=(${QUICKCHECK_CRAMS:-
    /mnt/8tb_hdd/flash_fraq_crams/HCC1395_SEQC2/CRAM/HCC1395_WGS_NV_3/HCC1395_WGS_NV_3.sorted.cram
    /mnt/8tb_hdd/flash_fraq_crams/HCC1395_SEQC2/CRAM/HCC1395BL_WGS_NV_3/HCC1395BL_WGS_NV_3.sorted.cram
})
# Extra --configfile overlays (space-separated), layered on top of the
# Snakefile's own configfile directives via Snakemake's recursive
# update_config -- see workflow/config/admixture_design_reshaped.yaml and
# workflow/config/admixture_reshaped.yaml for the reshaped run's overlay.
# shellcheck disable=SC2206
EXTRA_CONFIGFILES=(${EXTRA_CONFIGFILES:-})

mkdir -p logs
LOG="logs/admixture_mixtures_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1
echo "admixture_mixtures_local_driver started $(date -Is); log ${LOG}"

# --- pre-flight: source data is real, months-of-retrieval production data --
# never delete/regenerate it; this driver only ever reads it and only ever
# writes under ADMIX_OUTPUT_ROOT.

[ -s "${PLAN_PATH}" ] || { echo "ABORT: no admixture_plan.csv at ${PLAN_PATH} -- run the admixture_design checkpoint first." >&2; exit 1; }

echo "Checking the quickcheck source CRAMs/BAMs (${#QUICKCHECK_CRAMS[@]})..."
. "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate flash_fraq
for cram in "${QUICKCHECK_CRAMS[@]}"
do
    # A source may be a CRAM (.crai) or a BAM (.bai), e.g. the reshaped tumor BAM.
    [ -s "${cram}" ] && { [ -s "${cram}.crai" ] || [ -s "${cram}.bai" ]; } || { echo "ABORT: missing ${cram} or its .crai/.bai index" >&2; exit 1; }
    quickcheck_out="$(samtools quickcheck -v "${cram}" 2>&1)"
    [ -z "${quickcheck_out}" ] || { echo "ABORT: samtools quickcheck failed on ${cram}: ${quickcheck_out}" >&2; exit 1; }
    echo "  OK ${cram}"
done

echo "Checking every plan row's source path resolves..."
python3 - "${PLAN_PATH}" <<'PYEOF'
import csv, sys
plan_path = sys.argv[1]
import os
with open(plan_path) as f:
    rows = list(csv.DictReader(f))
missing = [
    (row["sample_name"], key, row[key])
    for row in rows
    for key in ("tumor_bam_path", "healthy_bam_path")
    if not os.path.exists(row[key])
]
print(f"  {len(rows)} plan rows, {len(missing)} missing source paths")
if missing:
    for m in missing[:20]:
        print("  MISSING", m)
    sys.exit(1)
PYEOF

echo "Checking free space on /mnt/8tb_hdd..."
AVAIL_GB="$(df -BG /mnt/8tb_hdd | tail -1 | awk '{print $4}' | tr -d 'G')"
echo "  ${AVAIL_GB}G available (need at least ${MIN_FREE_GB}G; each mixture's temp scratch is single-digit GB at C_target=1.0x and 4 run concurrently)"
if [ "${AVAIL_GB}" -lt "${MIN_FREE_GB}" ]; then
    echo "ABORT: not enough free space on /mnt/8tb_hdd." >&2
    exit 1
fi

conda activate snakemake
# localfeat_coverage shells out to mosdepth, which lives in its own conda env.
export PATH="${HOME}/miniconda3/envs/mosdepth/bin:${PATH}"
export PYTHONPATH="${REPO_ROOT}/src"

echo
echo "Lint..."
snakemake --lint -s "${SNAKEFILE}"

WORKFLOW_WORKDIR="${REPO_ROOT}/snakemake_workdirs/admixture"
mkdir -p "${WORKFLOW_WORKDIR}"
cd "${WORKFLOW_WORKDIR}"

SNAKEMAKE_ARGS=(-p -s "${SNAKEFILE}"
                --cores "${SNAKEMAKE_CORES}"
                --resources "mem_mb=${SNAKEMAKE_MEM_MB}"
                --keep-going
                # software-env excluded: rules declare conda: envs this driver does not
                # activate, so a --use-conda run elsewhere would otherwise re-trigger them.
                --rerun-triggers mtime params input code
                # rebuilds jobs left incomplete by an interrupted run.
                --rerun-incomplete)
if [ "${#EXTRA_CONFIGFILES[@]}" -gt 0 ]; then
    SNAKEMAKE_ARGS+=(--configfile "${EXTRA_CONFIGFILES[@]}")
    # "--" ends the multi-value --configfile list so the target is not read as a configfile.
fi

echo
echo "Dry run against the real target..."
snakemake -n "${SNAKEMAKE_ARGS[@]}" -- "${TARGET}"

echo
echo "Launching the real run: ${#SNAKEMAKE_ARGS[@]} base args, target ${TARGET}"
for pass in $(seq 1 "${PASSES}"); do
    echo
    echo "######## pass ${pass} of ${PASSES} -- $(date -Is) ########"
    if snakemake "${SNAKEMAKE_ARGS[@]}" -- "${TARGET}"; then
        echo "pass ${pass}: admixture_mixtures_complete.flag produced; nothing left to retry"
        break
    fi
    echo "pass ${pass}: finished with failures (see per-sample logs under ${ADMIX_OUTPUT_ROOT}/logs/admixture_mixture/)" >&2
    if [ "${pass}" -lt "${PASSES}" ]; then
        echo "sleeping 5m before the next pass"
        sleep 300
    fi
done

echo
echo "admixture_mixtures_local_driver finished $(date -Is)"
