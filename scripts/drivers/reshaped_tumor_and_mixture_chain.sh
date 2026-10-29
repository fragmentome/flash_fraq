#!/usr/bin/env bash
# Reshaped-tumor rebuild chain (docs/todo/todo/reshaped-mixture-downstream-rerun.md,
# docs/todo/in-progress/reshaping-pilot-post-hoc-filter.md's "go: full rebuild from
# pooled NV_1 + NV_3" verdict).
#
# Waits for the live upstream chain (upstream_chain_20260928.sh, 495 unreshaped
# mixtures) to finish, then: pools+dedups HCC1395_WGS_NV_1 + HCC1395_WGS_NV_3,
# measures its native coverage, reshapes it toward the pooled kept-healthy
# Cristiano fragment-size distribution at MAPQ>=30 genome-wide, patches the
# reshaped tumor's real native_coverage/bam_path into
# workflow/config/admixture_design_reshaped.yaml, then runs the admixture
# design + mixture-build DAG (workflow/rules/admixture_design.smk +
# admixture.smk) against the reshaped tumor source via the two overlay
# configfiles (admixture_design_reshaped.yaml, admixture_reshaped.yaml) --
# output lands under a brand-new root, admixture_mixtures_reshaped/, never
# colliding with the original admixture_mixtures/.
#
# matched_pair is deliberately NOT rebuilt: fragment-size-reshaping.md's own
# "reshaping doesn't apply to this regime" finding (no cfDNA reference
# distribution to mismatch a sheared-cell-line-vs-sheared-cell-line pair
# against). matched_pair_replicates_per_fraction=0 in the design overlay makes
# the new plan real_background-only; this script instead symlinks the 11
# existing matched_pair mixture directories from the original build into the
# new root so a later feature-matrix run over admixture_mixtures_reshaped/
# still has both regimes.
#
# Idempotent/resumable like upstream_chain_20260928.sh: every stage checks for
# its own already-good output first and skips. Heavy transient intermediates
# (the pooled/dedup tumor BAM) are staged on /mnt/data (NVMe, 593G free at
# last check) via a single piped samtools pipeline -- no multi-hundred-GB
# intermediate files ever hit disk simultaneously, only the two ~90G source
# CRAMs (read-only, untouched) and the final dedup BAM. Only small,
# already-reshaped output (~1x coverage, not ~238x) is written to
# /mnt/8tb_hdd; the 484 final mixture BAMs (the unavoidable ~1.6TB) are the
# only heavy writes this script sends to that disk, and only in the last
# stage, after the free-space gate below has passed.
#
# NOT launched automatically -- see this repo's prepared-work report
# (delivered via SubagentHandback) for the exact command to start this by
# hand, e.g.:
#   nohup ./scripts/drivers/reshaped_tumor_and_mixture_chain.sh \
#       > logs/reshaped_tumor_and_mixture_chain_$(date +%Y%m%d_%H%M%S).log 2>&1 &
set -euo pipefail
# Fail fast: the 2026-09-30 run kept going after the pool step died and still
# logged "CHAIN DONE". Any failing command now stops the chain.
trap 'echo "[$(date -Is)] reshaped_tumor_and_mixture_chain CHAIN FAILED at line ${LINENO}: ${BASH_COMMAND}" >&2' ERR
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"

# ---------------------------------------------------------------------------
# Config (env-overridable; defaults match the decided design)
# ---------------------------------------------------------------------------
UPSTREAM_LOG="${UPSTREAM_LOG:-/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures/logs/chain/upstream_chain_20260928.log}"
UPSTREAM_POLL_SECONDS="${UPSTREAM_POLL_SECONDS:-300}"

NV1_CRAM="${NV1_CRAM:-/mnt/8tb_hdd/flash_fraq_crams/HCC1395_SEQC2/CRAM/HCC1395_WGS_NV_1/HCC1395_WGS_NV_1.sorted.cram}"
NV3_CRAM="${NV3_CRAM:-/mnt/8tb_hdd/flash_fraq_crams/HCC1395_SEQC2/CRAM/HCC1395_WGS_NV_3/HCC1395_WGS_NV_3.sorted.cram}"
REFERENCE_FASTA="${REFERENCE_FASTA:-/mnt/data/reference_genome/hg19/hg19.fa}"

# Heavy transient intermediates: NVMe staging, never /mnt/8tb_hdd.
STAGING_ROOT="${STAGING_ROOT:-/mnt/data/flash_fraq_scratch/reshaped_tumor}"
# CRAM, not BAM: a ~238x BAM (~450G+) would not fit next to its own sort
# temp files on the 660G NVMe.
POOLED_DEDUP_BAM="${STAGING_ROOT}/HCC1395_NV1_NV3.pooled.dedup.cram"
POOL_PARTS_DIR="${STAGING_ROOT}/pool_parts"
POOL_JOBS="${POOL_JOBS:-5}"
POOL_SORT_THREADS="${POOL_SORT_THREADS:-3}"
POOL_SORT_MEM="${POOL_SORT_MEM:-2G}"
# Space-separated contig list for a quick test run; empty = every @SQ contig.
POOL_CONTIGS="${POOL_CONTIGS:-}"
RESHAPE_OUTPUT_DIR="${STAGING_ROOT}/reshaped"
RESHAPE_SAMPLE_NAME="HCC1395_NV1_NV3_pooled"
TARGET_HISTOGRAM="${TARGET_HISTOGRAM:-/mnt/8tb_hdd/flash_fraq_derived/reshaping_pilot/recheck_clean_regions/target_histogram.pooled.csv}"
MIN_MAPQ="${MIN_MAPQ:-30}"
# Per-contig process-pool fan-out for run_fragment_reshaping.py's histogram
# and retention passes (fragment_reshaping.compute_insert_size_histogram /
# reshape_bam). 24 leaves headroom on the 36-core box for the pool/dedup
# samtools pipeline's own -@ threads and mosdepth, should they ever overlap.
RESHAPE_WORKERS="${RESHAPE_WORKERS:-24}"

# Final, small output on /mnt/8tb_hdd -- must match
# workflow/config/admixture_design_reshaped.yaml's placeholder bam_path
# exactly (this script patches native_coverage/bam_path there after this
# lands).
FINAL_ADMIX_ROOT="${FINAL_ADMIX_ROOT:-/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures_reshaped}"
FINAL_RESHAPED_TUMOR_DIR="${FINAL_ADMIX_ROOT}/reshaped_tumor"
FINAL_RESHAPED_TUMOR_BAM="${FINAL_RESHAPED_TUMOR_DIR}/${RESHAPE_SAMPLE_NAME}.reshaped.bam"

DESIGN_OVERLAY="${REPO_ROOT}/workflow/config/admixture_design_reshaped.yaml"
ADMIX_OVERLAY="${REPO_ROOT}/workflow/config/admixture_reshaped.yaml"
PLAN_PATH="${FINAL_ADMIX_ROOT%/}/../cohort_reshaped/admixture_plan.csv"
PLAN_PATH="$(python3 -c 'import sys,pathlib; print(pathlib.Path(sys.argv[1]).resolve())' "${PLAN_PATH}")"

ORIGINAL_ADMIX_ROOT="${ORIGINAL_ADMIX_ROOT:-/mnt/8tb_hdd/flash_fraq_derived/admixture_mixtures}"
ORIGINAL_PLAN_PATH="${ORIGINAL_PLAN_PATH:-/mnt/8tb_hdd/flash_fraq_derived/cohort/admixture_plan.csv}"
BED_CLEANUP_SCRIPT="${ORIGINAL_ADMIX_ROOT}/logs/cleanup_fragment_beds.sh"

THREADS="${THREADS:-16}"
# mosdepth lives in its own conda env, not flash_fraq (see
# workflow/rules/local_features.smk's localfeat_coverage rule, `conda:
# "mosdepth"`) -- invoked by full path so the rest of this script can stay in
# the flash_fraq/snakemake envs without repeated conda activate churn.
MOSDEPTH_BIN="${MOSDEPTH_BIN:-$HOME/miniconda3/envs/mosdepth/bin/mosdepth}"

# Disk-space plan (see this script's companion report for the derivation):
# 484 real_background mixtures * ~3.3G + reshaped-tumor final copy (~1x
# coverage, small) on /mnt/8tb_hdd; generous NVMe headroom for the pool/dedup
# pipeline's own external-sort temp shards.
NEED_HDD_GB="${NEED_HDD_GB:-1650}"
NEED_DATA_GB="${NEED_DATA_GB:-250}"

FLASHFRAQ_CONDA_ENV="${FLASHFRAQ_CONDA_ENV:-flash_fraq}"
SNAKEMAKE_CONDA_ENV="${SNAKEMAKE_CONDA_ENV:-snakemake}"

mkdir -p "${REPO_ROOT}/logs"
LOG="${REPO_ROOT}/logs/reshaped_tumor_and_mixture_chain_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG}") 2>&1
ts(){ date -Is; }
echo "[$(ts)] reshaped_tumor_and_mixture_chain started; log ${LOG}"

. "$HOME/miniconda3/etc/profile.d/conda.sh"

# ---------------------------------------------------------------------------
# Stage 0: wait for the live upstream chain (495-mixture ichorCNA run) to
# finish. Never touches /mnt/8tb_hdd beyond `tail`ing a small log file.
# ---------------------------------------------------------------------------
echo "[$(ts)] waiting for CHAIN DONE in ${UPSTREAM_LOG}"
while ! grep -q "CHAIN DONE" "${UPSTREAM_LOG}" 2>/dev/null; do
    echo "[$(ts)] not done yet ($(tail -1 "${UPSTREAM_LOG}" 2>/dev/null)); sleeping ${UPSTREAM_POLL_SECONDS}s"
    sleep "${UPSTREAM_POLL_SECONDS}"
done
echo "[$(ts)] CHAIN DONE observed in ${UPSTREAM_LOG}; proceeding"

# ---------------------------------------------------------------------------
# Stage 1: free-space gate. Proposes (never executes without an explicit env
# opt-in) the already-prepared fragment-bed cleanup as the first thing to
# free space with, since it is reproducible-intermediate cleanup on data this
# same live chain just finished with.
# ---------------------------------------------------------------------------
avail_gb() { df -BG "$1" | tail -1 | awk '{print $4}' | tr -d 'G'; }

echo "[$(ts)] disk check: /mnt/8tb_hdd avail=$(avail_gb /mnt/8tb_hdd)G need>=${NEED_HDD_GB}G; /mnt/data avail=$(avail_gb /mnt/data)G need>=${NEED_DATA_GB}G"

if [ "${RUN_BED_CLEANUP:-0}" = "1" ] && [ -x "${BED_CLEANUP_SCRIPT}" ]; then
    echo "[$(ts)] RUN_BED_CLEANUP=1: executing ${BED_CLEANUP_SCRIPT} --execute (frees ~1.03T of reproducible fragments/*.bed on the now-finished 495-mixture set)"
    conda activate "${FLASHFRAQ_CONDA_ENV}"
    "${BED_CLEANUP_SCRIPT}" --execute
else
    echo "[$(ts)] RUN_BED_CLEANUP not set (or cleanup script missing) -- dry-run only, not freeing anything:"
    [ -x "${BED_CLEANUP_SCRIPT}" ] && "${BED_CLEANUP_SCRIPT}" --dry-run || echo "  (no cleanup script at ${BED_CLEANUP_SCRIPT})"
fi

HDD_AVAIL="$(avail_gb /mnt/8tb_hdd)"
DATA_AVAIL="$(avail_gb /mnt/data)"
if [ "${HDD_AVAIL}" -lt "${NEED_HDD_GB}" ]; then
    cat >&2 <<EOF
[$(ts)] ABORT: /mnt/8tb_hdd has only ${HDD_AVAIL}G free, need >= ${NEED_HDD_GB}G for the
484 real_background reshaped mixtures (~3.3G each) plus the reshaped tumor's
final copy. Freeing options, largest first (propose-only; none of these are
executed automatically):
  1. Run ${BED_CLEANUP_SCRIPT} --execute (RUN_BED_CLEANUP=1 here) -- frees
     ~1.03T of already-dead fragments/*.bed intermediates on the 495
     just-finished unreshaped mixtures (fragment_midpoints.bed,
     fragments_ichor_1mb.bed [zero readers, ever], fragments_delfi_5mb.bed --
     all reproducible from each mixture's dedup.bam).
  2. Once the downstream CNV-only interim model has consumed the feature
     matrix built from the *original* 495-mixture set, the 484
     real_background mix_* dedup BAMs under ${ORIGINAL_ADMIX_ROOT} become
     redundant with their reshaped replacements (~1.6T) -- archive or delete
     after that model's run is confirmed done and the reshaped rebuild has
     landed and been spot-checked. Keep the 11 matched_pair mix_* dirs
     (reused by this chain, not rebuilt). This is a real decision on
     production output -- propose only, do not act without the user.
EOF
    exit 1
fi
# The /mnt/data headroom is only for Stage 2 pooling; skip it once the pooled CRAM exists.
if [ ! -s "${POOLED_DEDUP_BAM}.crai" ] && [ "${DATA_AVAIL}" -lt "${NEED_DATA_GB}" ]; then
    echo "[$(ts)] ABORT: /mnt/data has only ${DATA_AVAIL}G free, need >= ${NEED_DATA_GB}G for the pool/dedup staging pipeline." >&2
    exit 1
fi
echo "[$(ts)] disk check passed"

conda activate "${FLASHFRAQ_CONDA_ENV}"
mkdir -p "${STAGING_ROOT}" "${RESHAPE_OUTPUT_DIR}" "${FINAL_RESHAPED_TUMOR_DIR}"

# ---------------------------------------------------------------------------
# Stage 2: pool + dedup NV_1 + NV_3 on /mnt/data, one contig at a time,
# POOL_JOBS contigs in parallel (largest first). Per contig:
#   merge -R contig -> namesort -> fixmate -m -> coordsort -> markdup -r -> CRAM
# then `samtools cat` the parts in @SQ header order (keeps coordinate order).
#
# Why per contig: the 2026-09-30 whole-genome pipe used `sort -u`, which also
# writes its temp shards uncompressed; at ~238x they filled the NVMe after
# 3.4 h ("failed writing ... Illegal seek"). Per-contig sorts with `-l 1`
# temp shards peak at a few tens of GB. Pairs split across contigs are
# orphaned by fixmate here, which is harmless: the reshaper keeps proper
# pairs only.
#
# `markdup -r` (removes duplicates), not merge_and_dedup's flag-only
# `markdup`: this file is a coverage/depth source for the reshaping math (the
# pilot's "C_dedup" figures are post-removal depths).
# ---------------------------------------------------------------------------
pool_contig() {
    local contig="$1" part="${POOL_PARTS_DIR}/$1.cram" tmp="${STAGING_ROOT}/tmp_pool/$1"
    if [ -s "${part}" ] && [ -e "${part}.ok" ]; then
        return 0
    fi
    mkdir -p "${tmp}"
    samtools merge -u -@ 2 --reference "${REFERENCE_FASTA}" -R "${contig}" - "${NV1_CRAM}" "${NV3_CRAM}" \
        | samtools sort -n -l 1 -m "${POOL_SORT_MEM}" -@ "${POOL_SORT_THREADS}" -T "${tmp}/ns" - \
        | samtools fixmate -m -u - - \
        | samtools sort -l 1 -m "${POOL_SORT_MEM}" -@ "${POOL_SORT_THREADS}" -T "${tmp}/cs" - \
        | samtools markdup -r -@ "${POOL_SORT_THREADS}" -T "${tmp}/md" \
            --reference "${REFERENCE_FASTA}" -O cram - "${part}.tmp"
    mv "${part}.tmp" "${part}"
    touch "${part}.ok"
    rm -rf "${tmp}"
    echo "[$(date -Is)] pooled ${contig}: $(du -h "${part}" | cut -f1)"
}
export -f pool_contig
export POOL_PARTS_DIR STAGING_ROOT REFERENCE_FASTA NV1_CRAM NV3_CRAM POOL_SORT_MEM POOL_SORT_THREADS

if [ -s "${POOLED_DEDUP_BAM}" ] && [ -s "${POOLED_DEDUP_BAM}.crai" ] && samtools quickcheck "${POOLED_DEDUP_BAM}" 2>/dev/null; then
    echo "[$(ts)] pooled/dedup tumor CRAM already present and passes quickcheck, skipping: ${POOLED_DEDUP_BAM}"
else
    mkdir -p "${POOL_PARTS_DIR}" "${STAGING_ROOT}/tmp_pool"
    if [ -n "${POOL_CONTIGS}" ]; then
        HEADER_CONTIGS="${POOL_CONTIGS}"
    else
        HEADER_CONTIGS="$(samtools view -H "${NV1_CRAM}" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2}' | tr '\n' ' ')"
    fi
    # Largest first, so chr1/chr2 don't start last and set the wall-clock.
    BY_SIZE="$(samtools view -H "${NV1_CRAM}" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); sub("LN:","",$3); print $3"\t"$2}' \
        | sort -k1,1nr | cut -f2 | grep -Fxf <(tr ' ' '\n' <<<"${HEADER_CONTIGS}" | grep -v '^$'))"
    echo "[$(ts)] pooling + dedup (removed) NV_1 + NV_3 per contig, ${POOL_JOBS} at a time -> ${POOL_PARTS_DIR}"
    printf '%s\n' ${BY_SIZE} | xargs -P "${POOL_JOBS}" -I{} bash -c 'pool_contig "$@"' _ {}
    PARTS=()
    for contig in ${HEADER_CONTIGS}; do
        [ -s "${POOL_PARTS_DIR}/${contig}.cram" ] && [ -e "${POOL_PARTS_DIR}/${contig}.cram.ok" ] \
            || { echo "[$(ts)] missing pool part for ${contig}" >&2; exit 1; }
        PARTS+=("${POOL_PARTS_DIR}/${contig}.cram")
    done
    echo "[$(ts)] concatenating ${#PARTS[@]} contig parts in header order"
    samtools cat -o "${POOLED_DEDUP_BAM}.tmp" "${PARTS[@]}"
    mv "${POOLED_DEDUP_BAM}.tmp" "${POOLED_DEDUP_BAM}"
    samtools index -@ "${THREADS}" "${POOLED_DEDUP_BAM}"
    samtools quickcheck "${POOLED_DEDUP_BAM}"
    # Parts are copies of the pooled CRAM's contents; drop them only now.
    rm -rf "${POOL_PARTS_DIR}" "${STAGING_ROOT}/tmp_pool"
    echo "[$(ts)] pool+dedup done: $(du -h "${POOLED_DEDUP_BAM}" | cut -f1)"
fi
# Smoke-test hook: STOP_AFTER_POOL=1 with POOL_CONTIGS/STAGING_ROOT set to
# scratch values exercises stage 2 alone.
if [ "${STOP_AFTER_POOL:-0}" = "1" ]; then
    echo "[$(ts)] STOP_AFTER_POOL=1, stopping after the pool stage"
    exit 0
fi

# ---------------------------------------------------------------------------
# Stage 3: measure the pooled/dedup source's native coverage (mosdepth,
# genome-wide mean, "total" row) -- purely diagnostic/sanity (expected
# ~238.6x per the pilot's tumor_depth/ check), not consumed downstream.
# ---------------------------------------------------------------------------
DEDUP_MOSDEPTH_PREFIX="${STAGING_ROOT}/HCC1395_NV1_NV3.pooled.dedup"
if [ -s "${DEDUP_MOSDEPTH_PREFIX}.mosdepth.summary.txt" ]; then
    echo "[$(ts)] pooled/dedup mosdepth summary already present, skipping"
else
    echo "[$(ts)] mosdepth on pooled/dedup source"
    "${MOSDEPTH_BIN}" -t "${THREADS}" -x --fasta "${REFERENCE_FASTA}" "${DEDUP_MOSDEPTH_PREFIX}" "${POOLED_DEDUP_BAM}"
fi
DEDUP_COVERAGE="$(awk '$1=="total"{print $4}' "${DEDUP_MOSDEPTH_PREFIX}.mosdepth.summary.txt")"
echo "[$(ts)] pooled/dedup native coverage: ${DEDUP_COVERAGE}x (pilot expected ~238.6x)"

# ---------------------------------------------------------------------------
# Stage 4: reshape genome-wide at MAPQ>=30 toward the pooled kept-healthy
# Cristiano target histogram (reused from the pilot, not recomputed --
# recheck_clean_regions/target_histogram.pooled.csv already *is* the
# genome-representative pooled target). One CLI call does histogram ->
# acceptance -> retention -> verification -> writes the reshaped BAM/index
# (scripts/run_fragment_reshaping.py, unmodified). Raises (nonzero exit) if
# the KS verification fails -- a real gate, not just advisory.
#
# NOTE: FragmentReshaper.run() is a per-read pysam Python loop, run twice
# full-genome over a ~238x-coverage BAM (source histogram pass + retention
# pass) plus once more over the much smaller output. At genome-wide depth
# this is orders of magnitude more reads than any pilot run to date (which
# only ever touched few-Mb regions). `--workers "${RESHAPE_WORKERS}"`
# (default 24) fans both full-BAM passes out per-contig across a process
# pool -- see fragment_reshaping.py's `compute_insert_size_histogram`/
# `reshape_bam` and this ticket's benchmark note in
# docs/todo/todo/reshaped-mixture-downstream-rerun.md for the measured
# per-worker throughput and the extrapolated wall-clock at 24 workers;
# pair accept/reject decisions are deterministic given (seed, qname, insert
# size), so the retained set does not depend on worker count.
# ---------------------------------------------------------------------------
if [ -s "${FINAL_RESHAPED_TUMOR_BAM}" ] && [ -s "${FINAL_RESHAPED_TUMOR_BAM}.bai" ] && samtools quickcheck "${FINAL_RESHAPED_TUMOR_BAM}" 2>/dev/null; then
    echo "[$(ts)] reshaped tumor BAM already present in its final location, skipping reshape: ${FINAL_RESHAPED_TUMOR_BAM}"
else
    STAGED_RESHAPED_BAM="${RESHAPE_OUTPUT_DIR}/${RESHAPE_SAMPLE_NAME}.reshaped.bam"
    if [ -s "${STAGED_RESHAPED_BAM}" ] && [ -s "${STAGED_RESHAPED_BAM}.bai" ] && samtools quickcheck "${STAGED_RESHAPED_BAM}" 2>/dev/null; then
        echo "[$(ts)] staged reshaped tumor BAM already present, skipping reshape step, will just copy: ${STAGED_RESHAPED_BAM}"
    else
        echo "[$(ts)] reshaping (MAPQ>=${MIN_MAPQ}) -> ${STAGED_RESHAPED_BAM}"
        python "${REPO_ROOT}/scripts/run_fragment_reshaping.py" \
            --source-bam "${POOLED_DEDUP_BAM}" \
            --reference-fasta "${REFERENCE_FASTA}" \
            --output-dir "${RESHAPE_OUTPUT_DIR}" \
            --target-histogram-path "${TARGET_HISTOGRAM}" \
            --min-mapq "${MIN_MAPQ}" \
            --sample-name "${RESHAPE_SAMPLE_NAME}" \
            --seed 42 \
            --workers "${RESHAPE_WORKERS}" \
            | tee "${RESHAPE_OUTPUT_DIR}/${RESHAPE_SAMPLE_NAME}.run_fragment_reshaping.json"
    fi
    echo "[$(ts)] copying reshaped tumor BAM to its final /mnt/8tb_hdd location: ${FINAL_RESHAPED_TUMOR_BAM}"
    cp "${STAGED_RESHAPED_BAM}" "${FINAL_RESHAPED_TUMOR_BAM}.tmp"
    mv "${FINAL_RESHAPED_TUMOR_BAM}.tmp" "${FINAL_RESHAPED_TUMOR_BAM}"
    cp "${STAGED_RESHAPED_BAM}.bai" "${FINAL_RESHAPED_TUMOR_BAM}.bai"
fi

# ---------------------------------------------------------------------------
# Stage 5: measure the *reshaped* BAM's own native coverage -- this is the
# real tumor.native_coverage the design's feasibility math needs (measured
# directly, not derived as C_dedup * ā, since the MAPQ filter also removes
# reads independently of the tlen-based acceptance draw).
# ---------------------------------------------------------------------------
RESHAPED_MOSDEPTH_PREFIX="${FINAL_RESHAPED_TUMOR_DIR}/${RESHAPE_SAMPLE_NAME}.reshaped"
if [ ! -s "${RESHAPED_MOSDEPTH_PREFIX}.mosdepth.summary.txt" ]; then
    echo "[$(ts)] mosdepth on reshaped tumor BAM"
    "${MOSDEPTH_BIN}" -t "${THREADS}" -x "${RESHAPED_MOSDEPTH_PREFIX}" "${FINAL_RESHAPED_TUMOR_BAM}"
fi
RESHAPED_COVERAGE="$(awk '$1=="total"{print $4}' "${RESHAPED_MOSDEPTH_PREFIX}.mosdepth.summary.txt")"
echo "[$(ts)] reshaped tumor native coverage: ${RESHAPED_COVERAGE}x (pilot's f_max~0.91 assumed ~0.907x-equivalent headroom at C_target=1.0x)"

# ---------------------------------------------------------------------------
# Stage 6: patch the design overlay's tumor block in place (native_coverage,
# bam_path) -- same file, same schema, only these two values change. Never
# touches admixture_design.yaml (the original/live config).
# ---------------------------------------------------------------------------
echo "[$(ts)] patching ${DESIGN_OVERLAY} tumor.native_coverage=${RESHAPED_COVERAGE}, tumor.bam_path=${FINAL_RESHAPED_TUMOR_BAM}"
# A targeted regex line-replace, not a yaml.safe_load/safe_dump round-trip:
# the overlay file's value is the whole point of it (see its own header
# comment) and a full re-dump would silently drop every comment in it.
# `count=1` on each substitution is what keeps this from also touching
# matched_pair_background's own (unrelated) native_coverage/bam_path keys --
# safe only because `tumor:` is written before `matched_pair_background:` in
# this file, so the first match of each key is always the tumor block's.
python3 - "${DESIGN_OVERLAY}" "${RESHAPED_COVERAGE}" "${FINAL_RESHAPED_TUMOR_BAM}" <<'PYEOF'
import re
import sys

overlay_path, native_coverage, bam_path = sys.argv[1], float(sys.argv[2]), sys.argv[3]
text = open(overlay_path).read()

text, n_cov = re.subn(
    r"^(\s*native_coverage:)\s*.*$",
    rf"\g<1> {native_coverage}",
    text,
    count=1,
    flags=re.MULTILINE,
)
text, n_bam = re.subn(
    r'^(\s*bam_path:)\s*.*$',
    rf'\g<1> "{bam_path}"',
    text,
    count=1,
    flags=re.MULTILINE,
)
if n_cov != 1 or n_bam != 1:
    raise SystemExit(
        f"expected exactly one native_coverage/bam_path line each in {overlay_path}, "
        f"matched {n_cov}/{n_bam} -- refusing to write a possibly-wrong patch"
    )
open(overlay_path, "w").write(text)
print(f"patched {overlay_path}: native_coverage={native_coverage}, bam_path={bam_path}")
PYEOF

# ---------------------------------------------------------------------------
# Stage 7: run the design checkpoint (native_coverage_manifest ->
# admixture_design), producing PLAN_PATH under cohort_reshaped/. Uses
# snakemake directly (not admixture_mixtures_local_driver.sh, which assumes
# the plan already exists) via the two overlay configfiles.
# ---------------------------------------------------------------------------
conda activate "${SNAKEMAKE_CONDA_ENV}"
export PYTHONPATH="${REPO_ROOT}/src"
DESIGN_WORKDIR="${REPO_ROOT}/snakemake_workdirs/admixture_design_reshaped"
mkdir -p "${DESIGN_WORKDIR}"
cd "${DESIGN_WORKDIR}"

# A plan older than the reshaped tumor BAM was designed against a different
# (or placeholder) tumor source -- the 2026-09-30 failed run left one -- so
# it is rebuilt with --allowed-rules + --forcerun on just the two design rules
# (staged CRAM mtimes would otherwise re-trigger all upstream local-feature jobs)
# rather than trusted.
if [ -s "${PLAN_PATH}" ] && [ "${PLAN_PATH}" -nt "${FINAL_RESHAPED_TUMOR_BAM}" ]; then
    echo "[$(ts)] admixture_plan.csv already present and newer than the reshaped tumor BAM, skipping design step: ${PLAN_PATH}"
else
    echo "[$(ts)] running admixture_design checkpoint -> ${PLAN_PATH}"
    snakemake -p -s "${REPO_ROOT}/workflow/Snakefile" \
        --cores 4 --resources mem_mb=8192 \
        --allowed-rules native_coverage_manifest admixture_design \
        --forcerun native_coverage_manifest admixture_design \
        --configfile "${DESIGN_OVERLAY}" "${ADMIX_OVERLAY}" \
        -- "${PLAN_PATH}"
fi
cd "${REPO_ROOT}"

# ---------------------------------------------------------------------------
# Stage 8: symlink the 11 existing matched_pair mixture directories from the
# original build into the new root, so a feature-matrix run over
# admixture_mixtures_reshaped/ has both regimes without rebuilding
# matched_pair (see header + matched_pair decision in the report). Reads the
# *original* plan (never modified) to find which mix_* rows are
# matched_pair; never touches ${ORIGINAL_ADMIX_ROOT} itself (symlinks only,
# read-only on the original).
# ---------------------------------------------------------------------------
echo "[$(ts)] symlinking matched_pair mixtures from the original build into ${FINAL_ADMIX_ROOT}"
python3 - "${ORIGINAL_PLAN_PATH}" "${ORIGINAL_ADMIX_ROOT}" "${FINAL_ADMIX_ROOT}" <<'PYEOF'
import csv
import sys
from pathlib import Path

original_plan_path, original_root, final_root = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
with open(original_plan_path) as f:
    rows = [r for r in csv.DictReader(f) if r["regime"] == "matched_pair"]

final_root.mkdir(parents=True, exist_ok=True)
linked, missing = 0, []
for row in rows:
    name = row["sample_name"]
    src = original_root / name
    dst = final_root / name
    if not src.is_dir():
        missing.append(str(src))
        continue
    if dst.exists() or dst.is_symlink():
        continue
    dst.symlink_to(src)
    linked += 1
print(f"matched_pair rows in original plan: {len(rows)}, newly symlinked: {linked}, missing sources: {len(missing)}")
for m in missing:
    print("  MISSING", m)
PYEOF

# ---------------------------------------------------------------------------
# Stage 9: build the real_background mixtures via the existing, now
# env-parameterised admixture_mixtures_local_driver.sh -- unchanged logic,
# new output_root/plan/quickcheck list/overlay configs only.
# ---------------------------------------------------------------------------
echo "[$(ts)] launching admixture_mixtures_local_driver.sh for the reshaped real_background build"
ADMIX_OUTPUT_ROOT="${FINAL_ADMIX_ROOT}" \
PLAN_PATH="${PLAN_PATH}" \
QUICKCHECK_CRAMS="${FINAL_RESHAPED_TUMOR_BAM} /mnt/8tb_hdd/flash_fraq_crams/HCC1395_SEQC2/CRAM/HCC1395BL_WGS_NV_3/HCC1395BL_WGS_NV_3.sorted.cram" \
EXTRA_CONFIGFILES="${DESIGN_OVERLAY} ${ADMIX_OVERLAY}" \
    "${REPO_ROOT}/scripts/drivers/admixture_mixtures_local_driver.sh"

echo "[$(ts)] reshaped_tumor_and_mixture_chain CHAIN DONE"
