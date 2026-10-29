#!/usr/bin/env bash
# Drive MAT-15 LOESS GC correction (src/flash_fraq/gc_correction.py, via
# scripts/run_local_gc_correction.py) over the 495 real admixture mixtures at
# /mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/ -- the real-data run
# gc-correction-loess.md's testing section has been waiting on since MAT-14's
# real fragment-extraction run completed 2026-09-21. This is what
# feature-matrix-assembler.md was actually blocked on.
#
# Not run through workflow/rules/local_features.smk's localfeat_gc_correct
# rule: that rule's input path is the cohort's per-dataset/sample_id derived
# layout. The mixtures sit in a flat mix_<name>/ layout on the HDD instead --
# same reasoning as admixture_mixture_fragments_driver.sh and
# admixture_mixture_ichor_driver.sh, this driver calls the tested script
# directly, one process per mixture, writing into the same mix_<name>/
# directory (gc_correction has no separate subdirectory of its own upstream
# either -- it writes alongside fragments' window_fragment_counts.tsv).
#
# No staging, no BAM: this rule operates purely on the small per-window TSV
# MAT-14 already wrote (pandas/statsmodels only), the one stage in the chain
# with no HDD-contention concern -- see run_local_gc_correction.py's
# docstring. Runs directly off the HDD.
#
# Concurrency: this is CPU-only pandas/statsmodels LOESS over a small table,
# not IO- or memory-bound like the other two mixture drivers -- 16 workers
# against 36 cores.
#
# Idempotent: skips any mixture whose arm_zscores.csv already exists, so a
# later pass only retries what failed or was interrupted.
#
# Fragment-bed cleanup: MAT-14's fragment_midpoints.bed/fragments_delfi_5mb.bed/
# fragments_ichor_1mb.bed (~2.5GB/mixture, ~1.1TB across all 495) are pure
# intermediates -- reproducible from the mixture BAM, and nothing past this
# rule reads any of them (fragments_ichor_1mb.bed has no reader at all; see
# flash_fraq.fragment_extraction.FRAGMENT_BED_INTERMEDIATES for the full
# audit). Once this driver's own gc-correction succeeds for a mixture, the
# three beds are dead weight, so this driver removes them by default --
# same default the Snakemake local_features.smk pipeline already has via
# temp(). Set KEEP_FRAGMENT_BEDS=1 to keep them (e.g. while debugging a
# fragment-extraction/GC-correction discrepancy).
#
# Usage:
#   nohup ./scripts/drivers/admixture_mixture_gc_correct_driver.sh \
#       > logs/admixture_mixture_gc_correct_YYYYMMDD_HHMMSS.log 2>&1 &

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"

# MIXTURES_ROOT is env-overridable (default unchanged) so this same,
# already-tested driver can also GC-correct the reshaped mixture set
# (admixture_mixtures_reshaped/) -- see
# scripts/drivers/reshaped_downstream_chain.sh.
MIXTURES_ROOT="${MIXTURES_ROOT:-/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures}"

JOBS="${JOBS:-16}"
# Coverage-bias mask (flash_fraq.genome_mask): drop 5 Mb windows within this many bp of
# an hg19 centromere or chromosome end before the LOESS fit. run_local_gc_correction.py
# defaults to the same value (genome_mask.DEFAULT_MASK_MARGIN_BP) when the flag is
# omitted; passed explicitly here only so this driver's behavior doesn't silently follow
# a future change to that default. 0 disables masking.
MASK_MARGIN_BP="${MASK_MARGIN_BP:-5000000}"

# See "Fragment-bed cleanup" above. 0 (default): remove the three fragment
# .bed intermediates after a successful gc-correction. 1: keep them.
KEEP_FRAGMENT_BEDS="${KEEP_FRAGMENT_BEDS:-0}"

mkdir -p logs
LOG="logs/admixture_mixture_gc_correct_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1
echo "admixture_mixture_gc_correct_driver started $(date -Is); log ${LOG}"

[ -d "${MIXTURES_ROOT}" ] || { echo "ABORT: missing ${MIXTURES_ROOT} -- run admixture_mixtures_local_driver.sh first." >&2; exit 1; }

. "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate flash_fraq

mapfile -t MIX_DIRS < <(find "${MIXTURES_ROOT}" -mindepth 1 -maxdepth 1 -type d -name 'mix_*' | sort)
echo "found ${#MIX_DIRS[@]} mixture directories"

run_one() {
    local mix_dir="$1"
    local sample_id
    sample_id="$(basename "${mix_dir}")"
    local window_counts="${mix_dir}/fragments/window_fragment_counts.tsv"
    if [ ! -s "${window_counts}" ]; then
        echo "SKIP ${sample_id}: no fragments/window_fragment_counts.tsv found"
        return 0
    fi
    local out_dir="${mix_dir}/fragments"
    local arm_zscores="${out_dir}/arm_zscores.csv"
    if [ -s "${arm_zscores}" ]; then
        echo "DONE ${sample_id}: already has arm_zscores.csv"
        return 0
    fi
    mkdir -p "${REPO_ROOT}/logs/admixture_mixture_gc_correct"
    local sample_log="${REPO_ROOT}/logs/admixture_mixture_gc_correct/${sample_id}.log"

    local cleanup_flag=()
    if [ "${KEEP_FRAGMENT_BEDS}" != "1" ]; then
        cleanup_flag=(--cleanup-fragment-beds)
    fi

    if PYTHONPATH="${REPO_ROOT}/src" python "${REPO_ROOT}/scripts/run_local_gc_correction.py" \
        --input-windows-path "${window_counts}" \
        --output-windows-path "${out_dir}/window_fragment_counts.gc_corrected.tsv" \
        --window-ratios-path "${out_dir}/window_ratios.csv" \
        --arm-zscore-path "${arm_zscores}" \
        --mask-margin-bp "${MASK_MARGIN_BP}" \
        "${cleanup_flag[@]}" \
        > "${sample_log}" 2>&1
    then
        echo "OK   ${sample_id}"
    else
        echo "FAIL ${sample_id}: see ${sample_log}"
    fi
}
export -f run_one
export REPO_ROOT MASK_MARGIN_BP KEEP_FRAGMENT_BEDS

printf '%s\n' "${MIX_DIRS[@]}" | xargs -I{} -P "${JOBS}" bash -c 'run_one "$@"' _ {}

total="${#MIX_DIRS[@]}"
done_count=$(find "${MIXTURES_ROOT}" -mindepth 3 -maxdepth 3 -name arm_zscores.csv -path '*/fragments/*' | wc -l)
echo "pass complete: ${done_count}/${total} mixtures have arm_zscores.csv"
