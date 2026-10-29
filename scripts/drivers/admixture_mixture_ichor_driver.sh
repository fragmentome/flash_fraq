#!/usr/bin/env bash
# Drive MAT-16 ichorCNA/HMMcopy binning (src/flash_fraq/cnv_binning.py, via
# scripts/run_local_cnv_binning.py) over the 495 real admixture mixtures at
# /mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/ -- the real-data run
# ichorcna-hmmcopy-wrapper.md's testing section has been waiting on since
# MAT-13's real run completed 2026-09-19. Modeled directly on
# admixture_mixture_fragments_driver.sh, which already did this for fragment
# extraction (2026-09-21).
#
# Not run through workflow/rules/local_features.smk's localfeat_ichor rule:
# that rule's inputs come from discover_local_crams over
# workflow/config/local_features.yaml's `datasets` (the cohort CRAM-mirror
# layout) and its cram_to_bam rule decodes a temp BAM from a staged CRAM. The
# mixtures are already BAMs (admixture_mixture's dedup step never produces a
# CRAM) sitting directly on the HDD in a flat mix_<name>/ layout -- this
# driver calls the tested script directly instead, one process per mixture,
# writing into an ichor/ subdirectory alongside each mixture's other outputs,
# mirroring the cohort pipeline's <dataset>/<sample_id>/ichor/ layout one
# level up.
#
# GC/mappability WIGs, centromere table, and ichorCNA script/libdir/binary
# paths are the same fixed values workflow/config/local_features.yaml pins
# for the cohort run (§8b.2's env split and the 0.3.2-vs-0.5.1 libdir
# mismatch), hardcoded here rather than re-parsed from the YAML since this
# driver runs outside Snakemake entirely.
#
# SSD staging: readCounter is BamTools-based, single-pass over the whole BAM;
# reused from the fragments driver's design to avoid the same HDD
# seek-contention problem local_features.smk's staging exists to avoid (see
# that file's docstring, §8b.4). Each mixture's .dedup.bam(+.bai) is rsynced
# to STAGING_ROOT before processing and removed after, regardless of success
# -- the HDD copy under MIXTURES_ROOT is never touched.
#
# Concurrency: 6 workers. runIchorCNA.R's EM fit is single-threaded but
# memory-heavy (16384 MB budgeted per sample in local_features.smk's
# localfeat_ichor rule) -- 6 workers is set by that memory budget
# (6*16GB=96GB against 103GB available as of 2026-09-24), not by CPU (36
# cores, single-threaded jobs) or HDD contention.
#
# Idempotent: skips any mixture whose ichor_estimates.csv already exists, so
# a later pass only retries what failed or was interrupted.
#
# Usage:
#   nohup ./scripts/drivers/admixture_mixture_ichor_driver.sh \
#       > logs/admixture_mixture_ichor_YYYYMMDD_HHMMSS.log 2>&1 &

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"

# MIXTURES_ROOT is env-overridable (default unchanged) so this same,
# already-tested driver can also run ichorCNA on the reshaped mixture set
# (admixture_mixtures_reshaped/) -- see
# scripts/drivers/reshaped_downstream_chain.sh.
MIXTURES_ROOT="${MIXTURES_ROOT:-/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures}"
REFERENCE_FASTA="/mnt/8tb_hdd/references/genome/hg19/hg19.fa"
STAGING_ROOT="${STAGING_ROOT:-/mnt/data/flash_fraq_stage/admixture_mixture_ichor}"

ICHORCNA_EXTDATA="/mnt/8tb_hdd/data/ichorCNA/inst/extdata"
GC_WIG="${ICHORCNA_EXTDATA}/gc_hg19_1000kb.wig"
MAP_WIG="${ICHORCNA_EXTDATA}/map_hg19_1000kb.wig"
CENTROMERE="${ICHORCNA_EXTDATA}/GRCh37.p13_centromere_UCSC-gapTable.txt"
ICHORCNA_SCRIPT="/mnt/8tb_hdd/data/ichorCNA/scripts/runIchorCNA.R"
ICHORCNA_LIBDIR="/mnt/8tb_hdd/data/ichorCNA"
NORMAL_PANEL="/mnt/8tb_hdd/data/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"
READCOUNTER_BIN="/home/matvii/miniconda3/envs/TSO500ichorCNA/bin/readCounter"
RSCRIPT_BIN="/home/matvii/miniconda3/envs/problem2/bin/Rscript"

CHROMOSOMES="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"
BIN_SIZE_BP=1000000
GENOME_BUILD="hg19"
GENOME_STYLE="UCSC"
ICHOR_PLOIDY="c(2,3)"
ICHOR_MAX_CN=5
ICHOR_TXN_E=0.9999
ICHOR_TXN_STRENGTH=10000

# Finding 8 of docs/results-corruption-audit.md: the mixtures are 1x ULP-WGS, and the
# previous default normal/scStates prior here was ichorCNA's *high-TF* recommendation
# (normal=c(0.5,...,0.9) plus subclonal states), which read flat/non-monotone TF at
# f<=0.05 (0.060, 0.048, 0.050, 0.065, 0.054 at f=0, 0.01, 0.02, 0.03, 0.05). ichorCNA's
# own low-TF-ULP-WGS recommendation is normal=c(0.95,0.99,0.995,0.999) with no
# subclonal states -- that is now the default. ICHOR_PROFILE=high_tf keeps the old
# values selectable, in case a future high-purity cohort needs them.
ICHOR_PROFILE="${ICHOR_PROFILE:-low_tf}"
case "${ICHOR_PROFILE}" in
    low_tf)
        ICHOR_NORMAL="c(0.95,0.99,0.995,0.999)"
        ICHOR_SC_STATES="c()"
        ICHOR_ESTIMATE_SC_PREVALENCE="FALSE"
        ;;
    high_tf)
        ICHOR_NORMAL="c(0.5,0.6,0.7,0.8,0.9)"
        ICHOR_SC_STATES="c(1,3)"
        ICHOR_ESTIMATE_SC_PREVALENCE="TRUE"
        ;;
    *)
        echo "ABORT: unknown ICHOR_PROFILE '${ICHOR_PROFILE}' (expected low_tf or high_tf)" >&2
        exit 1
        ;;
esac

JOBS="${JOBS:-6}"

mkdir -p logs
LOG="logs/admixture_mixture_ichor_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1
echo "admixture_mixture_ichor_driver started $(date -Is); log ${LOG}"

[ -s "${REFERENCE_FASTA}" ] || { echo "ABORT: missing ${REFERENCE_FASTA}" >&2; exit 1; }
[ -s "${GC_WIG}" ] || { echo "ABORT: missing ${GC_WIG}" >&2; exit 1; }
[ -s "${MAP_WIG}" ] || { echo "ABORT: missing ${MAP_WIG}" >&2; exit 1; }
[ -s "${CENTROMERE}" ] || { echo "ABORT: missing ${CENTROMERE}" >&2; exit 1; }
[ -s "${ICHORCNA_SCRIPT}" ] || { echo "ABORT: missing ${ICHORCNA_SCRIPT}" >&2; exit 1; }
[ -s "${NORMAL_PANEL}" ] || { echo "ABORT: missing ${NORMAL_PANEL}" >&2; exit 1; }
[ -x "${READCOUNTER_BIN}" ] || { echo "ABORT: missing ${READCOUNTER_BIN}" >&2; exit 1; }
[ -x "${RSCRIPT_BIN}" ] || { echo "ABORT: missing ${RSCRIPT_BIN}" >&2; exit 1; }
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
    local out_dir="${mix_dir}/ichor"
    if [ -s "${out_dir}/ichor_estimates.csv" ]; then
        echo "DONE ${sample_id}: already has ichor_estimates.csv"
        return 0
    fi
    mkdir -p "${REPO_ROOT}/logs/admixture_mixture_ichor"
    local sample_log="${REPO_ROOT}/logs/admixture_mixture_ichor/${sample_id}.log"

    local stage_dir="${STAGING_ROOT}/${sample_id}"
    local staged_bam="${stage_dir}/$(basename "${bam}")"
    mkdir -p "${stage_dir}"
    cleanup() { rm -rf "${stage_dir}"; }
    trap cleanup RETURN

    if ! rsync -a "${bam}" "${bam}.bai" "${stage_dir}/" > "${sample_log}" 2>&1; then
        echo "FAIL ${sample_id}: SSD staging copy failed, see ${sample_log}"
        return 0
    fi

    if PYTHONPATH="${REPO_ROOT}/src" python "${REPO_ROOT}/scripts/run_local_cnv_binning.py" \
        --alignment-path "${staged_bam}" \
        --reference-fasta "${REFERENCE_FASTA}" \
        --output-dir "${out_dir}" \
        --sample-id "${sample_id}" \
        --gc-wig "${GC_WIG}" \
        --map-wig "${MAP_WIG}" \
        --centromere "${CENTROMERE}" \
        --ichorcna-script "${ICHORCNA_SCRIPT}" \
        --ichorcna-libdir "${ICHORCNA_LIBDIR}" \
        --readcounter-bin "${READCOUNTER_BIN}" \
        --rscript-bin "${RSCRIPT_BIN}" \
        --bin-size-bp "${BIN_SIZE_BP}" \
        --genome-build "${GENOME_BUILD}" \
        --genome-style "${GENOME_STYLE}" \
        --chromosomes "${CHROMOSOMES}" \
        --normal-panel "${NORMAL_PANEL}" \
        --ploidy "${ICHOR_PLOIDY}" \
        --normal "${ICHOR_NORMAL}" \
        --sc-states "${ICHOR_SC_STATES}" \
        --estimate-sc-prevalence "${ICHOR_ESTIMATE_SC_PREVALENCE}" \
        --max-cn "${ICHOR_MAX_CN}" \
        --txn-e "${ICHOR_TXN_E}" \
        --txn-strength "${ICHOR_TXN_STRENGTH}" \
        >> "${sample_log}" 2>&1
    then
        echo "OK   ${sample_id}"
    else
        echo "FAIL ${sample_id}: see ${sample_log}"
    fi
}
export -f run_one
export REPO_ROOT REFERENCE_FASTA STAGING_ROOT GC_WIG MAP_WIG CENTROMERE
export ICHORCNA_SCRIPT ICHORCNA_LIBDIR NORMAL_PANEL READCOUNTER_BIN RSCRIPT_BIN
export CHROMOSOMES BIN_SIZE_BP GENOME_BUILD GENOME_STYLE
export ICHOR_PLOIDY ICHOR_NORMAL ICHOR_SC_STATES ICHOR_ESTIMATE_SC_PREVALENCE
export ICHOR_MAX_CN ICHOR_TXN_E ICHOR_TXN_STRENGTH

printf '%s\n' "${MIX_DIRS[@]}" | xargs -I{} -P "${JOBS}" bash -c 'run_one "$@"' _ {}

total="${#MIX_DIRS[@]}"
done_count=$(find "${MIXTURES_ROOT}" -mindepth 3 -maxdepth 3 -name ichor_estimates.csv -path '*/ichor/*' | wc -l)
echo "pass complete: ${done_count}/${total} mixtures have ichor_estimates.csv"
