#!/usr/bin/env bash
# Drive MAT-14 fragment extraction (src/flash_fraq/fragment_extraction.py, via
# scripts/run_local_fragment_extraction.py) over the 495 real admixture
# mixtures at /mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/ -- the "real
# merged BAM" run fragment-extraction.md's testing section has been waiting on
# since MAT-13's real run completed 2026-09-19.
#
# Not run through workflow/rules/local_features.smk's localfeat_fragments
# rule: that rule's inputs are CRAMs staged HDD->SSD via discover_local_crams
# over workflow/config/local_features.yaml's `datasets`, which is the cohort
# CRAM mirror layout (dataset/sample_id.sorted.cram). The mixtures are BAMs
# (admixture_mixture already produces a dedup .bam, no CRAM step) sitting
# directly on the HDD in a flat mix_<name>/ layout with no staging step --
# registering them as a `dataset` there would require re-plumbing the CRAM
# discovery/staging machinery for files that are neither CRAMs nor need
# staging. This driver calls the same tested script directly instead, one
# process per mixture, writing into a fragments/ subdirectory alongside each
# mixture's other outputs -- the same per-sample layout the cohort pipeline
# uses one level up.
#
# Window BEDs are shared with the cohort pipeline (genome-window definitions
# depend only on the hg19 .fai and primary-contig list, not on dataset) --
# reused from /mnt/8tb_hdd/flash_fraq_derived/windows/ rather than
# regenerated.
#
# SSD staging: measured live 2026-09-21 -- reading the mixture BAMs directly
# off the HDD with 4 concurrent workers ran at ~0.5 samples/min (>15h for all
# 495), well under the ~1 sample/min a single worker got in isolation. Same
# seek-contention problem local_features.smk's SSD-staging design exists to
# avoid (see that file's docstring, §8b.4). Each mixture's .dedup.bam(+.bai)
# is rsynced to STAGING_ROOT before processing and removed after, regardless
# of success -- the HDD copy under MIXTURES_ROOT is never touched.
#
# Concurrency: 8 workers. Each worker's HDD time is just its ~16s rsync copy
# (3GB dedup.bam); the ~2min pysam/bedtools extraction that follows is
# SSD+CPU-bound and doesn't touch the HDD, so at any instant only 1-2 of the
# 8 workers are typically mid-copy -- well short of the 3-4 concurrent
# streams §8b.4 found saturates a single HDD spindle. 36 cores comfortably
# cover 8 CPU-bound extractions. SSD headroom (603GB free as of 2026-09-21
# after the HCC1395 FASTQ cleanup) is not the constraint either.
#
# Idempotent: skips any mixture whose fragment_stats.csv already exists, so a
# later pass only retries what failed or was interrupted.
#
# Usage:
#   nohup ./scripts/drivers/admixture_mixture_fragments_driver.sh \
#       > logs/admixture_mixture_fragments_YYYYMMDD_HHMMSS.log 2>&1 &

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"

# MIXTURES_ROOT is env-overridable (default unchanged) so this same,
# already-tested driver can also run fragment extraction on the reshaped
# mixture set (admixture_mixtures_reshaped/) -- see
# scripts/drivers/reshaped_downstream_chain.sh.
MIXTURES_ROOT="${MIXTURES_ROOT:-/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures}"
WINDOWS_DIR="/mnt/8tb_hdd/flash_fraq_derived/windows"
DELFI_BED="${WINDOWS_DIR}/delfi_windows_5000000.bed"
ICHOR_BED="${WINDOWS_DIR}/ichor_bins_1000000.bed"
REFERENCE_FASTA="/mnt/8tb_hdd/references/genome/hg19/hg19.fa"
STAGING_ROOT="${STAGING_ROOT:-/mnt/data/flash_fraq_stage/admixture_mixture_fragments}"

JOBS="${JOBS:-8}"

mkdir -p logs
LOG="logs/admixture_mixture_fragments_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1
echo "admixture_mixture_fragments_driver started $(date -Is); log ${LOG}"

[ -s "${DELFI_BED}" ] || { echo "ABORT: missing ${DELFI_BED}" >&2; exit 1; }
[ -s "${ICHOR_BED}" ] || { echo "ABORT: missing ${ICHOR_BED}" >&2; exit 1; }
[ -s "${REFERENCE_FASTA}" ] || { echo "ABORT: missing ${REFERENCE_FASTA}" >&2; exit 1; }
[ -d "${MIXTURES_ROOT}" ] || { echo "ABORT: missing ${MIXTURES_ROOT} -- run admixture_mixtures_local_driver.sh first." >&2; exit 1; }

. "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate flash_fraq

mkdir -p "${STAGING_ROOT}"

mapfile -t MIX_DIRS < <(find "${MIXTURES_ROOT}" -mindepth 1 -maxdepth 1 -type d -name 'mix_*' | sort)
echo "found ${#MIX_DIRS[@]} mixture directories"

run_one() {
    local mix_dir="$1"
    local sample_id
    sample_id="$(basename "${mix_dir}")"
    local bam
    bam="$(find "${mix_dir}" -maxdepth 1 -name '*.dedup.bam' | head -1)"
    if [ -z "${bam}" ]; then
        echo "SKIP ${sample_id}: no .dedup.bam found"
        return 0
    fi
    local out_dir="${mix_dir}/fragments"
    if [ -s "${out_dir}/fragment_stats.csv" ]; then
        echo "DONE ${sample_id}: already has fragment_stats.csv"
        return 0
    fi
    mkdir -p "${REPO_ROOT}/logs/admixture_mixture_fragments"
    local sample_log="${REPO_ROOT}/logs/admixture_mixture_fragments/${sample_id}.log"

    local stage_dir="${STAGING_ROOT}/${sample_id}"
    local staged_bam="${stage_dir}/$(basename "${bam}")"
    mkdir -p "${stage_dir}"
    cleanup() { rm -rf "${stage_dir}"; }
    trap cleanup RETURN

    if ! rsync -a "${bam}" "${bam}.bai" "${stage_dir}/" > "${sample_log}" 2>&1; then
        echo "FAIL ${sample_id}: SSD staging copy failed, see ${sample_log}"
        return 0
    fi

    if PYTHONPATH="${REPO_ROOT}/src" python "${REPO_ROOT}/scripts/run_local_fragment_extraction.py" \
        --input-alignment "${staged_bam}" \
        --output-dir "${out_dir}" \
        --delfi-windows-bed "${DELFI_BED}" \
        --ichor-bins-bed "${ICHOR_BED}" \
        --reference-fasta "${REFERENCE_FASTA}" \
        --sample-id "${sample_id}" \
        >> "${sample_log}" 2>&1
    then
        echo "OK   ${sample_id}"
    else
        echo "FAIL ${sample_id}: see ${sample_log}"
    fi
}
export -f run_one
export REPO_ROOT DELFI_BED ICHOR_BED REFERENCE_FASTA STAGING_ROOT

printf '%s\n' "${MIX_DIRS[@]}" | xargs -I{} -P "${JOBS}" bash -c 'run_one "$@"' _ {}

total="${#MIX_DIRS[@]}"
done_count=$(find "${MIXTURES_ROOT}" -mindepth 3 -maxdepth 3 -name fragment_stats.csv -path '*/fragments/*' | wc -l)
echo "pass complete: ${done_count}/${total} mixtures have fragment_stats.csv"
