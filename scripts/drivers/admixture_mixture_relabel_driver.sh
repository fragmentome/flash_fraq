#!/usr/bin/env bash
# Re-run scripts/relabel_admixture_empirical_tf.py over all 495 real admixture
# mixtures to fix the empirical_tumor_fraction label bug in
# AdmixtureBuilder.compute_empirical_tumor_fraction: the merged/dedup BAM keeps
# the tumor and healthy sources' distinct @RG SM tags, so bcftools mpileup
# emitted one sample column per source and the old code read only the first
# (fields[9]), silently dropping the other source's reads at every site. Fixed
# in src/flash_fraq/admixture.py to sum AD across every sample column.
#
# No re-downsampling or re-merging needed: the dedup BAM on disk already
# reflects the correct mixture, only the label computed from it was wrong.
# This driver re-runs just the pileup + label step against the existing BAM
# and rewrites sample_metadata.csv + the per-label .txt files in place.
#
# Idempotent as of 2026-09-28: originally this always recomputed (the output
# always existed, it was just wrong, so "skip if output exists" didn't apply).
# Now that the fixed estimator also writes a `.sites.sites.tsv.gz` per-site
# dump as a side effect (see admixture.py's _write_empirical_tumor_fraction_sites),
# that dump's presence is a reliable signal that a mixture was already
# relabeled with the *current* code -- skip it on a later pass instead of
# redoing ~10 CPU-minutes of mpileup work. Delete the dump (or the whole
# mixture's sample_metadata.csv) to force a relabel.
#
# SSD staging added 2026-09-28: sda measured ~78% util at ~360 r/s / ~34MB/s
# during the first relabel pass -- seek-bound (bcftools mpileup -R over a
# sparse SNV site list forces many small indexed BAM reads scattered across
# the genome, unlike a single sequential scan). Same fix as
# admixture_mixture_fragments_driver.sh/admixture_mixture_ichor_driver.sh:
# rsync dedup.bam(+.bai) to STAGING_ROOT once (sequential HDD read), then run
# mpileup against the SSD copy, then clean up regardless of outcome.
#
# Usage:
#   nohup ./scripts/drivers/admixture_mixture_relabel_driver.sh \
#       > logs/admixture_mixture_relabel_YYYYMMDD_HHMMSS.log 2>&1 &

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"

# MIXTURES_ROOT is env-overridable (default unchanged) so this same,
# already-tested driver can also relabel the reshaped mixture set
# (admixture_mixtures_reshaped/) -- see
# scripts/drivers/reshaped_downstream_chain.sh.
MIXTURES_ROOT="${MIXTURES_ROOT:-/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures}"
REFERENCE_FASTA="/mnt/data/reference_genome/hg19/hg19.fa"
SNV_VAF_SITES="/mnt/data/reference_genome/hcc1395_seqc2/hcc1395_sSNV_hg19.sites.vcf.gz"
SNV_VAF_CLONAL="/mnt/data/reference_genome/hcc1395_seqc2/hcc1395_sSNV_hg19.clonal.vcf.gz"
STAGING_ROOT="${STAGING_ROOT:-/mnt/data/flash_fraq_stage/admixture_mixture_relabel}"

JOBS="${JOBS:-8}"

mkdir -p logs
LOG="logs/admixture_mixture_relabel_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1
echo "admixture_mixture_relabel_driver started $(date -Is); log ${LOG}"

[ -d "${MIXTURES_ROOT}" ] || { echo "ABORT: missing ${MIXTURES_ROOT}" >&2; exit 1; }

. "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate flash_fraq

mkdir -p "${STAGING_ROOT}"

mapfile -t MIX_DIRS < <(find "${MIXTURES_ROOT}" -mindepth 1 -maxdepth 1 -type d -name 'mix_*' | sort)
echo "found ${#MIX_DIRS[@]} mixture directories"

run_one() {
    local mix_dir="$1"
    local sample_id
    sample_id="$(basename "${mix_dir}")"
    local dedup_bam="${mix_dir}/${sample_id}.dedup.bam"
    if [ ! -s "${dedup_bam}" ]; then
        echo "SKIP ${sample_id}: no dedup BAM found"
        return 0
    fi
    local sites_dump="${mix_dir}/${sample_id}.empirical_tumor_fraction.sites.sites.tsv.gz"
    if [ -s "${sites_dump}" ]; then
        echo "DONE ${sample_id}: already relabeled with current estimator (has sites dump)"
        return 0
    fi
    mkdir -p "${REPO_ROOT}/logs/admixture_mixture_relabel"
    local sample_log="${REPO_ROOT}/logs/admixture_mixture_relabel/${sample_id}.log"

    local stage_dir="${STAGING_ROOT}/${sample_id}"
    local staged_bam="${stage_dir}/$(basename "${dedup_bam}")"
    mkdir -p "${stage_dir}"
    cleanup() { rm -rf "${stage_dir}"; }
    trap cleanup RETURN

    if ! rsync -a "${dedup_bam}" "${dedup_bam}.bai" "${stage_dir}/" > "${sample_log}" 2>&1; then
        echo "FAIL ${sample_id}: SSD staging copy failed, see ${sample_log}"
        return 0
    fi

    if PYTHONPATH="${REPO_ROOT}/src" python "${REPO_ROOT}/scripts/relabel_admixture_empirical_tf.py" \
        --mixture-dir "${mix_dir}" \
        --reference-fasta "${REFERENCE_FASTA}" \
        --snv-vaf-path "sites=${SNV_VAF_SITES}" \
        --snv-vaf-path "clonal=${SNV_VAF_CLONAL}" \
        --dedup-bam-path "${staged_bam}" \
        --threads 4 \
        >> "${sample_log}" 2>&1
    then
        echo "OK   ${sample_id}: $(tail -1 "${sample_log}")"
    else
        echo "FAIL ${sample_id}: see ${sample_log}"
    fi
}
export -f run_one
export REPO_ROOT REFERENCE_FASTA SNV_VAF_SITES SNV_VAF_CLONAL STAGING_ROOT

printf '%s\n' "${MIX_DIRS[@]}" | xargs -I{} -P "${JOBS}" bash -c 'run_one "$@"' _ {}

total="${#MIX_DIRS[@]}"
ok_count=$(grep -c '^OK' "${LOG}" || true)
fail_count=$(grep -c '^FAIL' "${LOG}" || true)
echo "pass complete: ${ok_count}/${total} OK, ${fail_count} FAIL"
