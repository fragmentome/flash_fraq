#!/usr/bin/env bash
# Take one arm of a SEQC2 tumor/normal pair, resolve the OTHER arm from the
# same sequencing site and library, and drive both to a verified CRAM on the
# 8 TB HDD -- then measure each arm's native coverage, which is the number
# admixture_design.yaml still carries as a placeholder.
#
# This is MAT-23's "ladder-unblocking set" deliverable end to end (fetch ->
# align -> verified CRAM -> coverage) for exactly one pair, which is all MAT-13
# needs to start building the fraction ladder.
#
# ## Why a pairing manifest exists at all
#
# MAT-23 requires that tumor and normal come from the SAME sequencing center --
# mixing sites reintroduces the batch confound that demoting Snyder to held-out
# validation was meant to remove. `center_name` in hcc1395_runs.csv cannot
# express that: it is `FDA` on all 118 rows, because FDA is the *submitting*
# center, not the sequencing site. Read as a site, it makes all 59 tumor rows
# look like a valid partner for any normal, and the confound comes back
# silently -- nothing downstream would flag it.
#
# The site lives only in ENA's `library_name`, which nothing else in this repo
# models. Its shape is <prep>_<site>_<arm>_<variant>:
#
#     WGS_IL_N_2   SRR7890858    56.3 Gbp      <- the aligned normal
#     WGS_IL_T_2   SRR7890863    67.8 Gbp      <- its pair: same site, same rep
#
# manifests/hcc1395/hcc1395_library_pairs.csv carries that decomposition for
# all 118 runs. SRP162370 holds four prep families -- WGS (73 runs), CHR (22),
# LBP (14, an input-mass titration) and SPP (9, pre-pooled genotype mixtures
# with no T/N arm at all) -- and this driver refuses to pair across them, since
# a CHR tumor against a WGS normal is a different confound wearing the same
# clothes.
#
# ## Dispatch is per arm, by what is already on disk
#
# Not a blind queue of two accessions. hcc1395_ena_local_driver.sh's
# clean_slate() deletes <run>.sorted.bam, so handing it a run whose alignment
# already finished would throw away hours of work and realign from scratch.
# Each arm therefore takes one of three paths:
#
#   CRAM + .crai on the HDD ......... nothing to do
#   sorted.bam, no CRAM ............. finish path: backdate the re-fetched
#                                     FASTQ below the BAM's mtime so
#                                     align_hcc1395 reads as up to date, then
#                                     run the ENA driver with
#                                     SKIP_CLEAN_SLATE=1 -- only cram_hcc1395
#                                     executes
#   nothing ......................... hand the whole run to the ENA driver
#
# Usage:
#   ./scripts/drivers/hcc1395_matched_pair_local_driver.sh                # SRR7890858
#   ./scripts/drivers/hcc1395_matched_pair_local_driver.sh SRR7890863
#   ./scripts/drivers/hcc1395_matched_pair_local_driver.sh --library WGS_IL_2
#   ./scripts/drivers/hcc1395_matched_pair_local_driver.sh --dry-run
#
# Env: SKIP_COVERAGE=1, MOSDEPTH_THREADS (4), plus everything
# hcc1395_ena_local_driver.sh reads (PREFETCH_AHEAD, MIN_FREE_GB,
# SNAKEMAKE_CORES, KEEP_SCRATCH, SKIP_HDD_PLACE, HDD_CRAM_ROOT ...), which is
# passed through untouched.
#
# Requires the same envs as the ENA driver (`snakemake`, `flash_fraq`) plus
# `mosdepth` for the coverage step.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"
ENA_DRIVER="${REPO_ROOT}/scripts/drivers/hcc1395_ena_local_driver.sh"
RUN_MANIFEST="${RUN_MANIFEST:-${REPO_ROOT}/manifests/hcc1395/hcc1395_runs.csv}"
# Overridable so the unmatched-arm guard below can be exercised: every T/N
# library in the real manifest has both arms, so nothing on disk can trigger it.
PAIR_MANIFEST="${PAIR_MANIFEST:-${REPO_ROOT}/manifests/hcc1395/hcc1395_library_pairs.csv}"

# Mirrors workflow/config/hcc1395_local.yaml, the same way the ENA driver does.
HCC1395_OUTPUT_ROOT="/mnt/data/flash_fraq_hcc1395"
REFERENCE_FASTA="/mnt/data/reference_genome/hg19/hg19.fa"
HDD_CRAM_ROOT="${HDD_CRAM_ROOT:-/mnt/8tb_hdd/flash_fraq_crams}"
HDD_DATASET="HCC1395_SEQC2"
BAM_ROOT="${HCC1395_OUTPUT_ROOT}/bam"
FASTQ_ROOT="${HCC1395_OUTPUT_ROOT}/fastq"
COVERAGE_ROOT="${HCC1395_OUTPUT_ROOT}/coverage"

MOSDEPTH_THREADS="${MOSDEPTH_THREADS:-4}"
ENA_FILEREPORT_URL="https://www.ebi.ac.uk/ena/portal/api/filereport"

usage() {
    cat <<'USAGE'
Usage: hcc1395_matched_pair_local_driver.sh [RUN | --library PREP_SITE_VARIANT] [--dry-run]

  RUN         either arm of the pair, e.g. SRR7890858; the other is resolved
  --library   name the library without its arm, e.g. WGS_IL_2 -> WGS_IL_N_2 + WGS_IL_T_2
  --dry-run   resolve and report the pair, then stop before any download

Env: SKIP_COVERAGE=1 MOSDEPTH_THREADS (4), plus every hcc1395_ena_local_driver.sh var
USAGE
    exit "${1:-0}"
}

# --- 1. parse args ----------------------------------------------------------

SEED_RUN=""
SEED_LIBRARY=""
DRY_RUN=0
while [ $# -gt 0 ]; do
    case "$1" in
        --library)  SEED_LIBRARY="$2"; shift 2 ;;
        --dry-run)  DRY_RUN=1; shift ;;
        -h|--help)  usage 0 ;;
        -*)         echo "ABORT: unknown flag $1" >&2; usage 1 >&2 ;;
        *)          [ -z "${SEED_RUN}" ] || { echo "ABORT: pass one accession, not two -- the second arm is resolved." >&2; exit 1; }
                    SEED_RUN="$1"; shift ;;
    esac
done
if [ -n "${SEED_RUN}" ] && [ -n "${SEED_LIBRARY}" ]; then
    echo "ABORT: --library and a run accession are two ways to say the same thing; pass one." >&2
    exit 1
fi
[ -n "${SEED_RUN}" ] || [ -n "${SEED_LIBRARY}" ] || SEED_RUN="${HCC1395_RUN:-SRR7890858}"

for f in "${RUN_MANIFEST}" "${PAIR_MANIFEST}" "${ENA_DRIVER}"; do
    [ -s "${f}" ] || { echo "ABORT: missing ${f}" >&2; exit 1; }
done

# hcc1395_runs.csv is CRLF (hcc1395_library_pairs.csv is not, but stripping a
# \r that is not there is free), so every awk pass over either drops it first.
STRIP_CR='{ sub(/\r$/, "") }'

in_run_manifest() {  # <run>: is this accession one of the 118 the chain knows?
    awk -F, -v run="$1" "${STRIP_CR} NR > 1 && \$1 == run { found = 1 } END { exit !found }" \
        "${RUN_MANIFEST}"
}

# Defined here rather than beside the dispatch code below because pair
# resolution already needs them: which run of a library to pick depends on
# what has been built for it.
hdd_run_dir() { echo "${HDD_CRAM_ROOT}/${HDD_DATASET}/CRAM/$1"; }

cram_done() {  # <run>
    local d; d="$(hdd_run_dir "$1")"
    [ -s "${d}/$1.sorted.cram" ] && [ -s "${d}/$1.sorted.cram.crai" ]
}

# --- 2. resolve the pair ----------------------------------------------------

# Fields of hcc1395_library_pairs.csv:
#   1 run_accession  2 library_name  3 prep  4 site  5 arm  6 variant
#   7 read_count     8 base_count
pair_field() {  # <run> <field-index>
    awk -F, -v run="$1" -v col="$2" "${STRIP_CR} NR > 1 && \$1 == run { print \$col; exit }" \
        "${PAIR_MANIFEST}"
}

# One run of a given library that the alignment chain actually knows about. A
# library maps to several runs -- WGS_IL_T_2 covers SRR7890863, SRR7910043,
# SRR7910044 and SRR7910045 (top-up lanes) -- so there is a choice to make.
#
# Work already on disk wins over depth, and by a long way: the deepest
# WGS_IL_N_2 run is SRR7909741 at 59.2 Gbp, but SRR7890858 (56.3 Gbp) was
# already aligned on 2026-09-03, and preferring the extra 2.9 Gbp would throw
# away a 32 GB sorted.bam and hours of BWA time for a difference the ladder
# downsamples away anyway. Only when no run of the library has been started
# does depth decide, since more depth is free headroom for the ladder.
best_run_for() {  # <prep> <site> <arm> <variant>
    local prep="$1" site="$2" arm="$3" variant="$4" run
    local started="" deepest=""
    while IFS=, read -r run _ p s a v _ _; do
        [ "${p}" = "${prep}" ] && [ "${s}" = "${site}" ] || continue
        [ "${a}" = "${arm}" ] && [ "${v}" = "${variant}" ] || continue
        in_run_manifest "${run}" || continue
        [ -n "${deepest}" ] || deepest="${run}"
        if [ -z "${started}" ] && { cram_done "${run}" || [ -s "${BAM_ROOT}/${run}.sorted.bam" ]; }; then
            started="${run}"
        fi
    done < <(awk -F, "${STRIP_CR} NR > 1" "${PAIR_MANIFEST}" | sort -t, -k8,8nr)
    if [ -n "${started}" ]; then echo "${started}"; return 0; fi
    [ -n "${deepest}" ] || return 1
    echo "${deepest}"
}

# ENA is the authority on library_name; the manifest is the offline fallback.
# Same live-then-committed order resolve_ena_pair() uses in the ENA driver, and
# for the same reason: the committed file must not silently go stale, but a
# flaky network must not stop the driver either.
library_name_for() {  # <run>
    local run="$1" live=""
    live="$(curl -fsS --retry 2 --max-time 30 \
        "${ENA_FILEREPORT_URL}?accession=${run}&result=read_run&fields=library_name&format=tsv" \
        2>/dev/null | awk -F'\t' 'NR == 2 { print $2 }' || true)"
    if [ -n "${live}" ]; then
        echo "${live}"
        return 0
    fi
    pair_field "${run}" 2
}

if [ -n "${SEED_LIBRARY}" ]; then
    # WGS_IL_2 -> prep WGS, site IL, variant 2. The arm is what we are adding.
    IFS=_ read -r SEED_PREP SEED_SITE SEED_VARIANT <<< "${SEED_LIBRARY}"
    if [ -z "${SEED_PREP}" ] || [ -z "${SEED_SITE}" ] || [ -z "${SEED_VARIANT}" ]; then
        echo "ABORT: --library wants <prep>_<site>_<variant>, e.g. WGS_IL_2 (got '${SEED_LIBRARY}')." >&2
        exit 1
    fi
else
    in_run_manifest "${SEED_RUN}" || {
        echo "ABORT: ${SEED_RUN} is not in ${RUN_MANIFEST}." >&2
        exit 1
    }
    SEED_LIB="$(library_name_for "${SEED_RUN}")"
    [ -n "${SEED_LIB}" ] || { echo "ABORT: no library_name for ${SEED_RUN}." >&2; exit 1; }
    SEED_PREP="$(pair_field "${SEED_RUN}" 3)"
    SEED_SITE="$(pair_field "${SEED_RUN}" 4)"
    SEED_ARM="$(pair_field "${SEED_RUN}" 5)"
    SEED_VARIANT="$(pair_field "${SEED_RUN}" 6)"
    if [ -z "${SEED_ARM}" ]; then
        echo "ABORT: ${SEED_RUN} (${SEED_LIB}) has no tumor/normal arm, so it has no pair." >&2
        echo "  SPP_* runs are pre-pooled genotype mixtures -- they are not a T/N pair." >&2
        exit 1
    fi
    # Guard against a manifest that disagrees with ENA about which library a
    # run belongs to: pairing off a stale row would pick the wrong site.
    if [ "${SEED_LIB}" != "${SEED_PREP}_${SEED_SITE}_${SEED_ARM}_${SEED_VARIANT}" ]; then
        echo "ABORT: ENA calls ${SEED_RUN} '${SEED_LIB}' but ${PAIR_MANIFEST} decomposes it as" >&2
        echo "  ${SEED_PREP}_${SEED_SITE}_${SEED_ARM}_${SEED_VARIANT}. Refresh the pairing manifest." >&2
        exit 1
    fi
fi

NORMAL_RUN="$(best_run_for "${SEED_PREP}" "${SEED_SITE}" N "${SEED_VARIANT}")" || NORMAL_RUN=""
TUMOR_RUN="$(best_run_for "${SEED_PREP}" "${SEED_SITE}" T "${SEED_VARIANT}")" || TUMOR_RUN=""

# An accession passed by hand is a decision, not a hint: only the arm the user
# did NOT name gets resolved. Without this, asking for SRR7890858 could hand
# back a different lane of the same library and quietly align the wrong run.
if [ -n "${SEED_RUN}" ]; then
    case "${SEED_ARM}" in
        N) NORMAL_RUN="${SEED_RUN}" ;;
        T) TUMOR_RUN="${SEED_RUN}" ;;
    esac
fi

# Never quietly substitute a different site or variant: an unmatched arm is the
# one failure this driver exists to make loud.
for arm in NORMAL TUMOR; do
    var="${arm}_RUN"
    if [ -z "${!var}" ]; then
        letter=N; [ "${arm}" = TUMOR ] && letter=T
        echo "ABORT: no ${arm,,} run for ${SEED_PREP}_${SEED_SITE}_${letter}_${SEED_VARIANT}" >&2
        echo "  in ${RUN_MANIFEST}. Pairing across sites or preps is exactly the confound" >&2
        echo "  this driver refuses to create -- pick a library with both arms." >&2
        exit 1
    fi
done

gbp() { awk -F, -v run="$1" "${STRIP_CR} NR > 1 && \$1 == run { printf \"%.1f\", \$8 / 1e9; exit }" "${PAIR_MANIFEST}"; }

echo
echo "=== matched pair: ${SEED_PREP}_${SEED_SITE}_*_${SEED_VARIANT} ==="
printf '  normal  %-12s %-14s %8s Gbp\n' "${NORMAL_RUN}" "$(pair_field "${NORMAL_RUN}" 2)" "$(gbp "${NORMAL_RUN}")"
printf '  tumor   %-12s %-14s %8s Gbp\n' "${TUMOR_RUN}" "$(pair_field "${TUMOR_RUN}" 2)" "$(gbp "${TUMOR_RUN}")"
echo "  final:  ${HDD_CRAM_ROOT}/${HDD_DATASET}/CRAM/<run>/"
echo

# --- 3. dispatch each arm ---------------------------------------------------

# Backdate both mates to one second before the sorted.bam. Snakemake compares
# mtimes, so a FASTQ newer than the BAM it produced makes align_hcc1395 look
# stale and triggers a rerun of the whole multi-hour alignment -- which is the
# trap that cost a realignment of SRR7890858 on 2026-09-03.
backdate_fastq() {  # <run>
    local run="$1" bam="${BAM_ROOT}/$1.sorted.bam" stamp mate
    stamp="$(date -d "@$(( $(stat -c%Y "${bam}") - 1 ))" '+%Y-%m-%d %H:%M:%S')"
    for mate in 1 2; do
        local path="${FASTQ_ROOT}/${run}/${run}_${mate}.fastq.gz"
        [ -f "${path}" ] || { echo "  ${run}: ${path##*/} is missing" >&2; return 1; }
        # gzip -t before trusting it: a still-resuming wget leaves a valid-
        # looking file that only fails at the end of the stream.
        gzip -t "${path}" 2>/dev/null || { echo "  ${run}: ${path##*/} is incomplete (gzip -t failed)" >&2; return 1; }
        touch -d "${stamp}" "${path}"
    done
    echo "  ${run}: backdated both mates to ${stamp} (sorted.bam is $(date -d "@$(stat -c%Y "${bam}")" '+%Y-%m-%d %H:%M:%S'))"
}

run_arm() {  # <run> <label>
    local run="$1" label="$2"
    echo "############ ${label}: ${run} ############"

    if cram_done "${run}"; then
        echo "  already at $(hdd_run_dir "${run}")/${run}.sorted.cram -- nothing to do"
        return 0
    fi

    if [ -s "${BAM_ROOT}/${run}.sorted.bam" ]; then
        echo "  sorted.bam present, no CRAM -- finish path (no realignment)"
        backdate_fastq "${run}" || {
            echo "  ${run}: cannot finish yet; let the FASTQ re-fetch complete, then re-run." >&2
            return 1
        }
        SKIP_CLEAN_SLATE=1 PREFETCH_AHEAD=0 "${ENA_DRIVER}" "${run}"
        return
    fi

    echo "  nothing on disk -- full fetch/align/cram through the ENA driver"
    "${ENA_DRIVER}" "${run}"
}

if [ "${DRY_RUN}" = "1" ]; then
    echo "--dry-run: stopping before any download."
    for pair in "normal:${NORMAL_RUN}" "tumor:${TUMOR_RUN}"; do
        label="${pair%%:*}"; run="${pair#*:}"
        if cram_done "${run}"; then      state="CRAM on HDD, nothing to do"
        elif [ -s "${BAM_ROOT}/${run}.sorted.bam" ]; then state="sorted.bam present -> finish path (cram only)"
        else                             state="nothing on disk -> full fetch/align/cram"
        fi
        printf '  %-7s %-12s %s\n' "${label}" "${run}" "${state}"
    done
    exit 0
fi

FAILED=()
# Normal first: it is the arm already part-built, so a failure there is cheap
# to see, and the tumor's multi-hour fetch does not run before we know it.
run_arm "${NORMAL_RUN}" normal || FAILED+=("${NORMAL_RUN}")
echo
run_arm "${TUMOR_RUN}" tumor || FAILED+=("${TUMOR_RUN}")

if [ ${#FAILED[@]} -gt 0 ]; then
    echo
    echo "FAILED: ${FAILED[*]} -- coverage skipped." >&2
    exit 1
fi

# --- 4. native coverage, the number admixture_design.yaml still fakes -------

if [ "${SKIP_COVERAGE:-0}" = "1" ]; then
    echo
    echo "SKIP_COVERAGE=1 -- both CRAMs are placed; coverage not measured."
    exit 0
fi

# Flags match rule localfeat_coverage (workflow/rules/local_features.smk), so
# these numbers sit on the same scale as the 547 rows already in
# native_coverage.csv rather than being a second, subtly different measurement.
measure_coverage() {  # <run> -> mean depth on stdout
    # Split, not one `local`: a single `local run="$1" prefix="...${run}"` has
    # its whole word list expanded BEFORE `local` assigns anything, so ${run}
    # resolves against the enclosing scope -- unbound under `set -u`, which is
    # exactly how this aborted after both CRAMs were already placed.
    local run="$1"
    local prefix="${COVERAGE_ROOT}/${run}" cram
    cram="$(hdd_run_dir "${run}")/${run}.sorted.cram"
    mkdir -p "${COVERAGE_ROOT}"
    if [ ! -s "${prefix}.mosdepth.summary.txt" ]; then
        conda run -n mosdepth mosdepth --threads "${MOSDEPTH_THREADS}" \
            --by 1000000 --no-per-base --fasta "${REFERENCE_FASTA}" \
            "${prefix}" "${cram}" >&2
    fi
    # The `total` row's mean, never `total_region` -- the same field and the
    # same reason as NativeCoverageHarvester.parse_mosdepth_summary().
    awk -F'\t' '$1 == "total" { print $4; found = 1; exit } END { exit !found }' \
        "${prefix}.mosdepth.summary.txt"
}

echo
echo "=== native coverage ==="
TUMOR_COV="$(measure_coverage "${TUMOR_RUN}")"
NORMAL_COV="$(measure_coverage "${NORMAL_RUN}")"

cat <<EOF

Paste into workflow/config/admixture_design.yaml, replacing the declared
placeholders (tumor 35.0 / matched_pair_background 34.0):

  tumor:
    native_coverage: ${TUMOR_COV}            # ${TUMOR_RUN} ($(pair_field "${TUMOR_RUN}" 2))
  matched_pair_background:
    native_coverage: ${NORMAL_COV}            # ${NORMAL_RUN} ($(pair_field "${NORMAL_RUN}" 2))

Note: bam_path in that config still points at sample-level CRAMs
(CRAM/HCC1395/HCC1395.sorted.cram), while this chain writes per-run
(CRAM/${TUMOR_RUN}/${TUMOR_RUN}.sorted.cram). Reconciling those, and adding
'${HDD_DATASET}' to local_features.yaml's datasets:, is MAT-23 Block C rollout
work -- not done here.
EOF
