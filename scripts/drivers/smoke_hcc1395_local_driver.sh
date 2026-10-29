#!/usr/bin/env bash
# Workstation twin of scripts/drivers/smoke_hcc1395_driver.sh: runs ONE SEQC2
# run accession through fetch_hcc1395_fastq -> align_hcc1395 -> cram_hcc1395
# on this machine, and times each stage
# (docs/todo/in-progress/hcc1395-local-smoke-benchmark.md).
#
# Differences from the cluster driver, and why:
#
#   * No ssh, no sbatch. Like smoke_local_features_driver.sh, this runs
#     Snakemake in place against the local checkout with
#     workflow/profiles/local, so an iteration costs no cluster round trip.
#   * Paths come from workflow/config/hcc1395_local.yaml, passed as a CLI
#     --configfile. Snakemake applies CLI configfiles after the Snakefile's
#     own `configfile:` directives and lets them win, so the cluster config
#     and every rule file stay untouched.
#   * The three rules are requested as three separate Snakemake invocations
#     rather than one .sorted.cram target. The DAG and the work are identical
#     -- each invocation just stops at the stage's own output -- but it yields
#     a per-stage wall clock, which a single invocation only buries in the
#     Snakemake log. This is a benchmark driver; the cluster one stays a
#     single target.
#
# Usage:
#   ./scripts/drivers/smoke_hcc1395_local_driver.sh
#
# Overridable via env var: HCC1395_RUN, MIN_FREE_GB, SKIP_CLEAN_SLATE=1.
#
# Requires the `snakemake` conda env (snakemake is not in `flash_fraq` --
# see docs/local-features-pipeline.md §8) and a `flash_fraq` env that actually
# satisfies environment.yml: sra-tools ships prefetch/fasterq-dump, and this
# box's env was missing it. The tool check below fails fast rather than
# letting a multi-hour run die at the first shell-out.

set -euo pipefail

# Same default accession as the cluster driver, so the two are comparable.
HCC1395_RUN="${HCC1395_RUN:-SRR7890858}"

cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"
SNAKEFILE="${REPO_ROOT}/workflow/Snakefile"
LOCAL_PROFILE="${REPO_ROOT}/workflow/profiles/local"
LOCAL_CONFIG="${REPO_ROOT}/workflow/config/hcc1395_local.yaml"

# Hardcoded to match workflow/config/hcc1395_local.yaml, the same way the
# cluster driver mirrors workflow/config/hcc1395.yaml: the target paths and
# the free-space check are built from them.
HCC1395_OUTPUT_ROOT="/mnt/data/flash_fraq_hcc1395"
REFERENCE_FASTA="/mnt/data/reference_genome/hg19/hg19.fa"
# fasterq-dump writes ~150G of uncompressed FASTQ for one SEQC2 run before
# gzip, on top of the .sra and an unsorted+sorted BAM pair. Check before
# starting rather than discovering it mid-run -- an exhausted fileset is what
# failed MAT-23's cluster run.
MIN_FREE_GB="${MIN_FREE_GB:-300}"

RUN_FASTQ_DIR="${HCC1395_OUTPUT_ROOT}/fastq/${HCC1395_RUN}"
FASTQ_TARGET="${RUN_FASTQ_DIR}/${HCC1395_RUN}_1.fastq.gz"
BAM_TARGET="${HCC1395_OUTPUT_ROOT}/bam/${HCC1395_RUN}.sorted.bam"
CRAM_TARGET="${HCC1395_OUTPUT_ROOT}/bam/${HCC1395_RUN}.sorted.cram"

. "$HOME/miniconda3/etc/profile.d/conda.sh"

echo "Checking the flash_fraq env has the tools the rules shell out to..."
missing_tools=()
for tool in prefetch fasterq-dump bwa samtools gzip; do
    conda run -n flash_fraq command -v "${tool}" > /dev/null 2>&1 || missing_tools+=("${tool}")
done
if [ ${#missing_tools[@]} -gt 0 ]; then
    echo "ABORT: missing from the flash_fraq env: ${missing_tools[*]}" >&2
    echo "Refresh it against environment.yml:  conda env update -f environment.yml --prune" >&2
    exit 1
fi

echo "Checking the local hg19 reference and its BWA index..."
missing_ref=()
for suffix in "" .fai .amb .ann .bwt .pac .sa; do
    [ -s "${REFERENCE_FASTA}${suffix}" ] || missing_ref+=("${REFERENCE_FASTA}${suffix}")
done
if [ ${#missing_ref[@]} -gt 0 ]; then
    echo "ABORT: missing or empty reference file(s): ${missing_ref[*]}" >&2
    echo "Fetch them first:  ./scripts/drivers/fetch_hg19_reference_local_driver.sh" >&2
    exit 1
fi

mkdir -p "${HCC1395_OUTPUT_ROOT}"
AVAIL_GB="$(df -BG "${HCC1395_OUTPUT_ROOT}" | tail -1 | awk '{print $4}' | tr -d 'G')"
echo "  ${AVAIL_GB}G available on ${HCC1395_OUTPUT_ROOT} (need at least ${MIN_FREE_GB}G)"
if [ "${AVAIL_GB}" -lt "${MIN_FREE_GB}" ]; then
    echo "ABORT: not enough free space for this run (~150G of uncompressed FASTQ before gzip)." >&2
    exit 1
fi

if [ "${SKIP_CLEAN_SLATE:-0}" != "1" ]; then
    echo
    echo "Preparing a clean slate for ${HCC1395_RUN}..."
    # Same reasoning as the cluster driver: a run interrupted during the gzip
    # step leaves <run>_2.fastq next to a fresh <run>.sra_2.fastq, and the
    # retriever's mate-resolution glob then matches both and aborts -- after
    # the multi-hour extraction. Abandoned fasterq.tmp.* scratch is dropped
    # for the same reason it was on isilon: it is dead weight on the disk.
    # The .sra is KEPT so a re-run's prefetch reports "found locally" and
    # skips re-downloading tens of GB.
    rm -f "${RUN_FASTQ_DIR}"/*_[12].fastq
    rm -rf "${RUN_FASTQ_DIR}"/fasterq.tmp.*
    rm -f "${BAM_TARGET}" "${BAM_TARGET}.bai"
    ls -la "${RUN_FASTQ_DIR}" 2>/dev/null || echo "  (no ${RUN_FASTQ_DIR} yet)"
fi

conda activate snakemake
export PYTHONPATH="${REPO_ROOT}/src"

# Own .snakemake -- lock, metadata, logs -- separate from the cluster
# workflows' and from local_features' (snakemake_workdirs/README.md), so this
# can run while another workflow holds its own directory. Every path handed
# to Snakemake above is absolute for that reason.
WORKFLOW_WORKDIR="${REPO_ROOT}/snakemake_workdirs/hcc1395_local"
mkdir -p "${WORKFLOW_WORKDIR}"
cd "${WORKFLOW_WORKDIR}"

# --configfile takes a LIST, so it must not be the last flag: with the target
# path directly after it, Snakemake swallows the target as a second config
# file and dies on "No such file or directory: <target>".
# --cores overrides the local profile's `cores: 16`, which exists to stop
# concurrent local_features samples from thrashing the disk. Nothing runs
# concurrently here -- it is one sample, three stages -- so the cap would just
# hand bwa mem 16 of this box's 36 cores and halve the number this benchmark
# is meant to produce. Must match the config's `threads:` (32), since
# Snakemake silently clamps a rule's threads to the available cores.
SNAKEMAKE_CORES="${SNAKEMAKE_CORES:-32}"
SNAKEMAKE_ARGS=(-p --use-conda -s "${SNAKEFILE}"
                --configfile "${LOCAL_CONFIG}"
                --workflow-profile "${LOCAL_PROFILE}"
                --cores "${SNAKEMAKE_CORES}")

echo
echo "Pre-flight checks..."
snakemake --lint -s "${SNAKEFILE}"
snakemake -n "${SNAKEMAKE_ARGS[@]}" "${CRAM_TARGET}"

STAGE_NAMES=(fetch align cram)
STAGE_TARGETS=("${FASTQ_TARGET}" "${BAM_TARGET}" "${CRAM_TARGET}")
STAGE_SECONDS=()

for i in "${!STAGE_NAMES[@]}"; do
    echo
    echo "=== ${STAGE_NAMES[$i]}: ${STAGE_TARGETS[$i]} ==="
    stage_start="${SECONDS}"
    snakemake "${SNAKEMAKE_ARGS[@]}" "${STAGE_TARGETS[$i]}"
    STAGE_SECONDS+=("$(( SECONDS - stage_start ))")
done

echo
echo "=== Wall clock, ${HCC1395_RUN} on $(hostname) ($(nproc) cores) ==="
total=0
for i in "${!STAGE_NAMES[@]}"; do
    secs="${STAGE_SECONDS[$i]}"
    total=$(( total + secs ))
    printf '  %-6s %6ds  (%dh%02dm)\n' "${STAGE_NAMES[$i]}" "${secs}" "$(( secs / 3600 ))" "$(( secs % 3600 / 60 ))"
done
printf '  %-6s %6ds  (%dh%02dm)\n' total "${total}" "$(( total / 3600 ))" "$(( total % 3600 / 60 ))"

echo
echo "Artifacts:"
ls -la "${CRAM_TARGET}" "${CRAM_TARGET}.crai" 2>&1
echo "(expect no ${BAM_TARGET} and no ${RUN_FASTQ_DIR}/${HCC1395_RUN}.sra left behind)"
