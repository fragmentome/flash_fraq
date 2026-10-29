#!/usr/bin/env bash
# Take one or more run accessions from manifests/hcc1395/, pull their FASTQ
# pairs straight from ENA, align each to a verified CRAM on this workstation,
# and land those CRAMs on the 8 TB HDD in the mirror layout local_features
# already consumes.
#
# Why this exists next to smoke_hcc1395_local_driver.sh, which runs the same
# three rules: that driver's fetch stage goes through HCC1395FastqRetriever,
# which only consults ENA's *submitted* files and falls back to SRA
# prefetch/fasterq-dump when a run has none. For SRP162370 submitted_ftp is
# empty on every run -- SEQC2 submitted BAMs -- so every run takes the slow
# path: a ~26 GB .sra staged, then ~150 GB of uncompressed FASTQ before gzip.
# ENA's *generated* fastq.gz pair, which manifests/hcc1395/hcc1395_ena_wget.sh
# already enumerates for all 118 runs, is ~35 GB and needs no round trip. This
# driver downloads that pair itself and lets Snakemake pick the chain up at
# align_hcc1395, so fetch_hcc1395_fastq never runs.
#
# Nothing under workflow/ or src/ changes for this: the download writes the
# exact paths fetch_hcc1395_fastq declares as its outputs, so Snakemake sees
# that rule as already satisfied. Fixing the retriever to consult fastq_ftp is
# worth doing for the cluster too, but that is a separate, tested change.
#
# ## Downloads run ahead of the alignment, not in front of it
#
# Measured against ENA on 2026-09-03: one connection sustains ~3-6 MB/s, four
# concurrent connections sustain ~2.4-3.2 MB/s EACH (~11 MB/s aggregate) on a
# 1 Gbit link that is nowhere near saturated. ENA throttles per connection, so
# concurrency is close to free. Two consequences, both implemented here:
#
#   * A run's R1 and R2 download in parallel, halving a pair's ~90 min.
#   * A background prefetcher works down the queue while the main loop aligns.
#     Alignment is CPU-bound on all 32 cores and uses no network; downloading
#     uses the network and almost no CPU. So every run after the first finds
#     its FASTQ already on disk and the fetch stage leaves the critical path
#     entirely.
#
# The prefetcher is bounded by two things, because the queue can be 118 runs x
# 35 GB and the NVMe holds ~590 GB: PREFETCH_AHEAD ready-but-unconsumed pairs,
# and a free-space floor that also has to leave room for the aligner's own
# ~120 GB of BAM scratch. It blocks rather than fills the disk.
#
# Usage:
#   ./scripts/drivers/hcc1395_ena_local_driver.sh                 # default run
#   ./scripts/drivers/hcc1395_ena_local_driver.sh SRR7890858 SRR7890859
#   ./scripts/drivers/hcc1395_ena_local_driver.sh --group hcc1395bl
#   ./scripts/drivers/hcc1395_ena_local_driver.sh --all
#
# Overridable via env var: PREFETCH_AHEAD (2; 0 disables prefetch, a large
# value downloads as far ahead as the disk budget allows), MIN_FREE_GB (250),
# MIN_FREE_HDD_GB (50), SNAKEMAKE_CORES (32), KEEP_SCRATCH=1, SKIP_HDD_PLACE=1,
# SKIP_CLEAN_SLATE=1 (keep an already-verified sorted.bam -- see clean_slate),
# HDD_CRAM_ROOT.
#
# Requires the `snakemake` conda env (snakemake is not in `flash_fraq` -- see
# docs/local-features-pipeline.md §8) and a `flash_fraq` env with bwa and
# samtools. No SRA toolkit is needed on this path.

set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/../.."
REPO_ROOT="${PWD}"
SNAKEFILE="${REPO_ROOT}/workflow/Snakefile"
LOCAL_PROFILE="${REPO_ROOT}/workflow/profiles/local"
LOCAL_CONFIG="${REPO_ROOT}/workflow/config/hcc1395_local.yaml"
RUN_MANIFEST="${REPO_ROOT}/manifests/hcc1395/hcc1395_runs.csv"
ENA_WGET_MANIFEST="${REPO_ROOT}/manifests/hcc1395/hcc1395_ena_wget.sh"

# Hardcoded to match workflow/config/hcc1395_local.yaml, the same way the other
# local driver mirrors it: the target paths and the free-space check are built
# from these. The NVMe is scratch only -- the FASTQ pair, the unsorted BAM and
# the sorted BAM all live and die there; only the CRAM reaches the HDD.
HCC1395_OUTPUT_ROOT="/mnt/data/flash_fraq_hcc1395"
REFERENCE_FASTA="/mnt/data/reference_genome/hg19/hg19.fa"
# <dataset>/CRAM/<sample>/<sample>.sorted.cram is what sync_crams_to_hdd.sh
# produces from the cluster archive and what local_features.yaml's `cram_root`
# reads, so a locally produced CRAM lands in the same shape as a mirrored one.
HDD_CRAM_ROOT="${HDD_CRAM_ROOT:-/mnt/8tb_hdd/flash_fraq_crams}"
HDD_DATASET="HCC1395_SEQC2"

# Floor kept free on the NVMe. One run needs ~35 GB of FASTQ plus an unsorted
# and a sorted BAM (~60 GB each for a 56 Gbp run), and the prefetcher must not
# eat the room the aligner is about to need.
MIN_FREE_GB="${MIN_FREE_GB:-250}"
MIN_FREE_HDD_GB="${MIN_FREE_HDD_GB:-50}"
PREFETCH_AHEAD="${PREFETCH_AHEAD:-2}"

ENA_FILEREPORT_URL="https://www.ebi.ac.uk/ena/portal/api/filereport"

FASTQ_ROOT="${HCC1395_OUTPUT_ROOT}/fastq"
BAM_ROOT="${HCC1395_OUTPUT_ROOT}/bam"
FETCH_LOG_DIR="${HCC1395_OUTPUT_ROOT}/logs/ena_fetch"
# Markers coordinating the background prefetcher with the align loop. Kept out
# of the fastq tree so nothing here can be mistaken for a rule output.
STATE_DIR="${HCC1395_OUTPUT_ROOT}/.ena_driver_state"

usage() {
    cat <<'USAGE'
Usage: hcc1395_ena_local_driver.sh [RUN...] [--all] [--group G] [--tumor-status S]

  RUN...            run accessions from manifests/hcc1395/hcc1395_runs.csv
  --all             every run in that manifest
  --group G         every run with source_group G (hcc1395 | hcc1395bl)
  --tumor-status S  every run with tumor_status S (tumor | healthy)

Env: PREFETCH_AHEAD (2, 0 disables) MIN_FREE_GB (250) MIN_FREE_HDD_GB (50)
     SNAKEMAKE_CORES (32) KEEP_SCRATCH=1 SKIP_HDD_PLACE=1 SKIP_CLEAN_SLATE=1
     HDD_CRAM_ROOT
USAGE
    exit "${1:-0}"
}

# --- 1. build the run queue from the committed manifest ---------------------

[ -s "${RUN_MANIFEST}" ] || { echo "ABORT: missing ${RUN_MANIFEST}" >&2; exit 1; }

# hcc1395_runs.csv is CRLF, so without the sub() below the last field
# (tumor_status) carries a trailing \r and `$6 == "tumor"` never matches.
# Assigning to $0 re-splits the fields, so it has to come before any $N test.
STRIP_CR='{ sub(/\r$/, "") }'

manifest_runs() {  # all accessions, in manifest order
    awk -F, "${STRIP_CR} NR > 1 && \$1 != \"\" { print \$1 }" "${RUN_MANIFEST}"
}
manifest_row() {  # <run> -> the manifest line, empty if unknown
    awk -F, -v run="$1" "${STRIP_CR} NR > 1 && \$1 == run { print; exit }" "${RUN_MANIFEST}"
}

QUEUE=()
selected=()
while [ $# -gt 0 ]; do
    case "$1" in
        --all)    mapfile -t selected < <(manifest_runs)
                  QUEUE+=("${selected[@]}"); shift ;;
        --group)  mapfile -t selected < <(awk -F, -v g="$2" "${STRIP_CR} NR > 1 && \$5 == g { print \$1 }" "${RUN_MANIFEST}")
                  [ ${#selected[@]} -gt 0 ] || { echo "ABORT: no runs with source_group '$2'." >&2; exit 1; }
                  QUEUE+=("${selected[@]}"); shift 2 ;;
        --tumor-status)
                  mapfile -t selected < <(awk -F, -v s="$2" "${STRIP_CR} NR > 1 && \$6 == s { print \$1 }" "${RUN_MANIFEST}")
                  [ ${#selected[@]} -gt 0 ] || { echo "ABORT: no runs with tumor_status '$2'." >&2; exit 1; }
                  QUEUE+=("${selected[@]}"); shift 2 ;;
        -h|--help) usage 0 ;;
        -*)       echo "ABORT: unknown flag $1" >&2; usage 1 >&2 ;;
        *)        QUEUE+=("$1"); shift ;;
    esac
done
[ ${#QUEUE[@]} -gt 0 ] || QUEUE=("${HCC1395_RUN:-SRR7890858}")

for run in "${QUEUE[@]}"; do
    if [ -z "$(manifest_row "${run}")" ]; then
        echo "ABORT: ${run} is not in ${RUN_MANIFEST}." >&2
        echo "Known accessions:" >&2
        awk -F, "${STRIP_CR} NR > 1 { print \"  \" \$1 \"  \" \$2 \"  \" \$6 }" "${RUN_MANIFEST}" >&2
        exit 1
    fi
done

hdd_run_dir() { echo "${HDD_CRAM_ROOT}/${HDD_DATASET}/CRAM/$1"; }
already_done() {  # <run>: CRAM + index already on the HDD
    local d; d="$(hdd_run_dir "$1")"
    [ -s "${d}/$1.sorted.cram" ] && [ -s "${d}/$1.sorted.cram.crai" ]
}

# Drop what is already on the HDD before anything is downloaded for it.
PENDING=()
for run in "${QUEUE[@]}"; do
    if already_done "${run}"; then
        echo "skip ${run}: already at $(hdd_run_dir "${run}")/${run}.sorted.cram"
    else
        PENDING+=("${run}")
    fi
done
if [ ${#PENDING[@]} -eq 0 ]; then
    echo "Nothing to do -- every requested run is already on the HDD."
    exit 0
fi

echo
echo "=== queue: ${#PENDING[@]} run(s) ==="
for run in "${PENDING[@]}"; do
    IFS=, read -r _ sample _ platform group status <<< "$(manifest_row "${run}")"
    printf '  %-12s %-10s %-10s %-8s %s\n' "${run}" "${sample}" "${group}" "${status}" "${platform}"
done
echo "  prefetch ahead: ${PREFETCH_AHEAD}   scratch: ${HCC1395_OUTPUT_ROOT}"
echo "  final:          ${HDD_CRAM_ROOT}/${HDD_DATASET}/CRAM/<run>/"
echo

# --- 2. pre-flight ----------------------------------------------------------

. "$HOME/miniconda3/etc/profile.d/conda.sh"

echo "Checking the flash_fraq env has the tools the rules shell out to..."
missing_tools=()
# No prefetch/fasterq-dump: this path never reaches the SRA toolkit.
for tool in bwa samtools; do
    conda run -n flash_fraq command -v "${tool}" > /dev/null 2>&1 || missing_tools+=("${tool}")
done
for tool in wget curl md5sum rsync flock; do
    command -v "${tool}" > /dev/null 2>&1 || missing_tools+=("${tool}")
done
if [ ${#missing_tools[@]} -gt 0 ]; then
    echo "ABORT: missing tool(s): ${missing_tools[*]}" >&2
    echo "For the conda ones, refresh the env:  conda env update -f environment.yml --prune" >&2
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

mkdir -p "${FASTQ_ROOT}" "${BAM_ROOT}" "${FETCH_LOG_DIR}" "${STATE_DIR}"

free_gb() { df -BG "$1" | tail -1 | awk '{print $4}' | tr -d 'G'; }

AVAIL_GB="$(free_gb "${HCC1395_OUTPUT_ROOT}")"
echo "  ${AVAIL_GB}G free on ${HCC1395_OUTPUT_ROOT} (need ${MIN_FREE_GB}G)"
if [ "${AVAIL_GB}" -lt "${MIN_FREE_GB}" ]; then
    echo "ABORT: not enough scratch space (~35G FASTQ + an unsorted/sorted BAM pair)." >&2
    exit 1
fi

if [ "${SKIP_HDD_PLACE:-0}" != "1" ]; then
    # Refuse to write into a bare mountpoint directory if the HDD is not
    # mounted -- that would silently fill the little space left on /, the same
    # failure sync_crams_to_hdd.sh guards against.
    HDD_MOUNT="$(echo "${HDD_CRAM_ROOT}" | cut -d/ -f1-3)"
    mountpoint -q "${HDD_MOUNT}" || { echo "ABORT: ${HDD_MOUNT} is not mounted." >&2; exit 1; }
    HDD_AVAIL_GB="$(free_gb "${HDD_MOUNT}")"
    echo "  ${HDD_AVAIL_GB}G free on ${HDD_MOUNT} (need ${MIN_FREE_HDD_GB}G)"
    if [ "${HDD_AVAIL_GB}" -lt "${MIN_FREE_HDD_GB}" ]; then
        echo "ABORT: not enough room on ${HDD_MOUNT} for the CRAM." >&2
        exit 1
    fi
fi

# One driver at a time: two would fight over the same Snakemake work directory
# (DAG build dies on the lock) and over the same scratch budget.
exec 9> "${STATE_DIR}/driver.lock"
if ! flock -n 9; then
    echo "ABORT: another hcc1395_ena_local_driver.sh is running (${STATE_DIR}/driver.lock)." >&2
    exit 1
fi

# --- 3. seed the run manifest the checkpoint would otherwise resolve --------

# align_hcc1395 takes <output_root>/manifest/hcc1395_runs.csv (the
# resolve_hcc1395_runs checkpoint's output) as an input, and common.smk's
# _run_metadata() reads source_group/tumor_status/platform out of it. Copying
# the committed manifest into place satisfies the checkpoint, so no online SRA
# query happens here.
SEEDED_MANIFEST="${HCC1395_OUTPUT_ROOT}/manifest/hcc1395_runs.csv"
mkdir -p "$(dirname "${SEEDED_MANIFEST}")"
if ! cmp -s "${RUN_MANIFEST}" "${SEEDED_MANIFEST}"; then
    cp "${RUN_MANIFEST}" "${SEEDED_MANIFEST}"
    echo "Seeded ${SEEDED_MANIFEST} from ${RUN_MANIFEST}."
fi

# --- 4. the ENA fetch side --------------------------------------------------

clean_slate() {  # <run>: drop debris from the SRA path this driver replaces
    local run="$1" dir="${FASTQ_ROOT}/$1"
    # A stalled prefetch leaves <run>.sra.prf/.sra.tmp (16 GB was sitting here
    # from 2026-08-31), and a run interrupted during gzip leaves an
    # uncompressed <run>_2.fastq that the retriever's mate glob would later
    # trip over. The fastq.gz pair is deliberately KEPT -- fetch resumes it.
    rm -f "${dir}/${run}".sra "${dir}/${run}".sra.prf "${dir}/${run}".sra.tmp
    rm -f "${dir}"/*_[12].fastq
    rm -rf "${dir}"/fasterq.tmp.* "${dir}/${run}".sra.tmp.*
    # The BAMs are dropped so a re-run realigns from scratch rather than
    # trusting a half-written sort. SKIP_CLEAN_SLATE=1 keeps them, which is the
    # one case where that is wrong-headed: a run whose sorted.bam is already
    # verified and whose FASTQ was merely re-fetched to satisfy the DAG needs
    # only cram_hcc1395, and deleting the BAM here would cost a realignment of
    # several hours. hcc1395_matched_pair_local_driver.sh's finish path sets it
    # after backdating the FASTQ mtimes; the dry-run guard below still refuses
    # to continue if fetch/resolve are somehow back in the DAG.
    if [ "${SKIP_CLEAN_SLATE:-0}" != "1" ]; then
        rm -f "${BAM_ROOT}/${run}.sorted.bam" "${BAM_ROOT}/${run}.sorted.bam.bai" \
              "${BAM_ROOT}/${run}.unsorted.bam"
    fi
}

resolve_ena_pair() {  # <run> -> "R1_URL<TAB>R2_URL<TAB>R1_MD5<TAB>R2_MD5"
    local run="$1" tsv fields ftp md5
    tsv="$(curl -fsS --retry 3 --max-time 60 \
        "${ENA_FILEREPORT_URL}?accession=${run}&result=read_run&fields=fastq_ftp,fastq_md5,fastq_bytes&format=tsv" \
        2>/dev/null || true)"
    fields="$(printf '%s\n' "${tsv}" | awk -F'\t' 'NR == 2 { print $2 "\t" $3 }')"
    ftp="$(printf '%s' "${fields}" | cut -f1)"
    md5="$(printf '%s' "${fields}" | cut -f2)"
    if [ -n "${ftp}" ]; then
        printf '%s\t%s\t%s\t%s\n' \
            "https://$(printf '%s' "${ftp}" | cut -d';' -f1)" \
            "https://$(printf '%s' "${ftp}" | cut -d';' -f2)" \
            "$(printf '%s' "${md5}" | cut -d';' -f1)" \
            "$(printf '%s' "${md5}" | cut -d';' -f2)"
        return 0
    fi
    # ENA unreachable or serving nothing for this run: fall back to the URLs
    # already committed in the wget manifest. No md5s come with those.
    printf '%s\t%s\t\t\n' \
        "$(grep -o "ftp://[^ ]*${run}_1\.fastq\.gz" "${ENA_WGET_MANIFEST}" | head -n1)" \
        "$(grep -o "ftp://[^ ]*${run}_2\.fastq\.gz" "${ENA_WGET_MANIFEST}" | head -n1)"
}

verify_fastq() {  # <path> <expected md5 or empty>
    local path="$1" expected="$2"
    # 1 MB floor: gzip -t already rejects a truncated download, but an empty
    # gzip stream is a valid 20-byte file, and treating one as "done" would
    # send wget -c resuming from byte 20 into a corrupt member.
    [ -f "${path}" ] && [ "$(stat -c%s "${path}")" -gt 1000000 ] || return 1
    if [ -n "${expected}" ]; then
        [ "$(md5sum "${path}" | cut -d' ' -f1)" = "${expected}" ] || return 1
    fi
    gzip -t "${path}"
}

fetch_fastq() {  # <url> <path> <expected md5 or empty> <wget log>
    local url="$1" path="$2" expected="$3" wlog="$4"
    if verify_fastq "${path}" "${expected}" 2>/dev/null; then
        echo "    ${path##*/}: already present and verified"
        return 0
    fi
    # -c resumes a partial file, and the URL's basename is already the name the
    # fetch_hcc1395_fastq rule declares, so -P lands it straight at its final
    # path -- no rename step that a killed download could leave half-applied.
    # -o, not the terminal: the two mates download concurrently and their dot
    # progress would otherwise interleave into noise.
    #
    # The retry budget is sized for scheduled network maintenance, not for a
    # flaky link: MedUni Graz IT announced short outages on 2026-09-07/08/10
    # between 20:00 and 24:00, and a fetch that exhausts its tries here is not
    # retried -- fetch_run_pair marks the run <run>.failed and the align loop
    # skips it for good. The old --tries=20 --waitretry=15 gave ~3 min of
    # linear backoff, which a "wenige Minuten" outage can outlast. 200 tries at
    # up to 60 s apart is ~2.9 h of budget. A genuinely dead URL still fails
    # fast: wget does not retry 404/403 no matter how large --tries is.
    # --retry-on-host-error added 2026-09-10, after SRR7909924_2 died on
    # "unable to resolve host address" at try 3 of 200 during the announced
    # 20:00-24:00 window. --tries governs connect and read errors only; wget
    # treats a DNS resolution failure as fatal and exits at once, so the large
    # budget above did not apply to it at all. This flag reclassifies host
    # errors as transient, putting them under --tries like every other failure.
    wget -c --tries=200 --waitretry=60 --retry-connrefused --retry-on-host-error \
         --timeout=60 \
         --progress=dot:giga \
         -o "${wlog}" -P "$(dirname "${path}")" "${url}" || {
        echo "    ${path##*/}: wget failed, see ${wlog}" >&2
        return 1
    }
    if ! verify_fastq "${path}" "${expected}"; then
        echo "    ${path##*/}: FAILED md5/gzip verification" >&2
        return 1
    fi
    echo "    ${path##*/}: verified"
}

fetch_run_pair() {  # <run>: both mates in parallel; 0 on success
    local run="$1" dir="${FASTQ_ROOT}/$1"
    local r1="${dir}/${run}_1.fastq.gz" r2="${dir}/${run}_2.fastq.gz"
    mkdir -p "${dir}"
    clean_slate "${run}"

    local urls r1_url r2_url r1_md5 r2_md5
    urls="$(resolve_ena_pair "${run}")"
    IFS=$'\t' read -r r1_url r2_url r1_md5 r2_md5 <<< "${urls}"
    if [ -z "${r1_url}" ] || [ -z "${r2_url}" ]; then
        echo "  ${run}: no ENA FASTQ pair URL" >&2
        return 1
    fi
    [ -n "${r1_md5}" ] || echo "  ${run}: no md5s from ENA -- gzip integrity only"

    # ENA throttles per connection (~3-6 MB/s), so the two mates together move
    # about twice what one does. Their exit codes are collected separately;
    # `wait` on a failed job must not take the whole prefetcher down.
    local rc=0 p1 p2
    fetch_fastq "${r1_url}" "${r1}" "${r1_md5}" "${FETCH_LOG_DIR}/${run}_1.wget.log" & p1=$!
    fetch_fastq "${r2_url}" "${r2}" "${r2_md5}" "${FETCH_LOG_DIR}/${run}_2.wget.log" & p2=$!
    wait "${p1}" || rc=1
    wait "${p2}" || rc=1
    [ "${rc}" -eq 0 ] || return 1

    # Same shape as RetrievalRecord.to_dict() in
    # src/flash_fraq/hcc1395_retrieval.py, so a run retrieved this way stays as
    # auditable as one the retriever fetched.
    cat > "${dir}/${run}_retrieval.json" <<JSON
{
  "run_accession": "${run}",
  "source": "ena_generated_fastq",
  "r1_path": "${r1}",
  "r2_path": "${r2}",
  "source_urls": ["${r1_url}", "${r2_url}"]
}
JSON
}

# --- 5. the background prefetcher -------------------------------------------

# Marker protocol, all in STATE_DIR:
#   <run>.fetched   written by the prefetcher; contains the fetch seconds
#   <run>.failed    written by the prefetcher; the align loop skips that run
#   <run>.consumed  written by the align loop when it starts on a run, so the
#                   prefetcher knows that pair no longer counts against the
#                   ready-ahead budget (the aligner deletes the FASTQ itself)
rm -f "${STATE_DIR}"/*.fetched "${STATE_DIR}"/*.failed "${STATE_DIR}"/*.consumed

ready_ahead() {  # pairs downloaded but not yet picked up by the align loop
    local n=0 f run
    for f in "${STATE_DIR}"/*.fetched; do
        [ -e "${f}" ] || continue
        run="$(basename "${f}" .fetched)"
        [ -e "${STATE_DIR}/${run}.consumed" ] || n=$((n + 1))
    done
    echo "${n}"
}

prefetch_worker() {
    local run start
    for run in "${PENDING[@]}"; do
        # Two brakes, both re-checked every 30 s: how far ahead of the aligner
        # this is allowed to get, and how much room is left on the NVMe once
        # the aligner's own BAM scratch is accounted for.
        while [ "$(ready_ahead)" -ge "${PREFETCH_AHEAD}" ] || \
              [ "$(free_gb "${HCC1395_OUTPUT_ROOT}")" -lt "${MIN_FREE_GB}" ]; do
            sleep 30
        done
        echo "  [prefetch] ${run}: downloading R1+R2"
        start="${SECONDS}"
        if fetch_run_pair "${run}"; then
            echo "$(( SECONDS - start ))" > "${STATE_DIR}/${run}.fetched"
            echo "  [prefetch] ${run}: ready in $(( SECONDS - start ))s"
        else
            touch "${STATE_DIR}/${run}.failed"
            echo "  [prefetch] ${run}: FAILED" >&2
        fi
    done
}

kill_tree() {  # <pid>: the worker's wgets are grandchildren, so pkill -P alone
    local pid="$1" child        # would leave them downloading after a Ctrl-C
    for child in $(pgrep -P "${pid}" 2>/dev/null); do kill_tree "${child}"; done
    kill "${pid}" 2>/dev/null || true
}

PREFETCH_PID=""
if [ "${PREFETCH_AHEAD}" -gt 0 ]; then
    prefetch_worker &
    PREFETCH_PID=$!
    # Do not leave a detached downloader running if the driver is interrupted.
    trap '[ -n "${PREFETCH_PID}" ] && kill_tree "${PREFETCH_PID}"; true' EXIT
fi

# --- 6. align -> cram -> HDD, one run at a time -----------------------------

conda activate snakemake
export PYTHONPATH="${REPO_ROOT}/src"

# Own .snakemake -- lock, metadata, logs -- separate from hcc1395_local's
# (snakemake_workdirs/README.md), so this can run while another workflow holds
# its own directory. Every path handed to Snakemake is absolute for that reason.
WORKFLOW_WORKDIR="${REPO_ROOT}/snakemake_workdirs/hcc1395_ena_local"
mkdir -p "${WORKFLOW_WORKDIR}"
cd "${WORKFLOW_WORKDIR}"

# --configfile takes a LIST, so it must not be the last flag: with the target
# path directly after it, Snakemake swallows the target as a second config file
# and dies on "No such file or directory: <target>".
SNAKEMAKE_CORES="${SNAKEMAKE_CORES:-32}"
SNAKEMAKE_ARGS=(-p --use-conda -s "${SNAKEFILE}"
                --configfile "${LOCAL_CONFIG}"
                --workflow-profile "${LOCAL_PROFILE}"
                --cores "${SNAKEMAKE_CORES}")

echo
echo "Pre-flight checks..."
snakemake --lint -s "${SNAKEFILE}"

SUMMARY=()
FAILED=()

for run in "${PENDING[@]}"; do
    bam_target="${BAM_ROOT}/${run}.sorted.bam"
    cram_target="${BAM_ROOT}/${run}.sorted.cram"
    run_dir="$(hdd_run_dir "${run}")"

    echo
    echo "############ ${run} ############"

    if [ "${PREFETCH_AHEAD}" -gt 0 ]; then
        waited="${SECONDS}"
        while [ ! -e "${STATE_DIR}/${run}.fetched" ] && [ ! -e "${STATE_DIR}/${run}.failed" ]; do
            # The prefetcher dying (OOM, kill) would otherwise hang this loop
            # forever waiting on a marker nobody is going to write.
            if ! kill -0 "${PREFETCH_PID}" 2>/dev/null && \
               [ ! -e "${STATE_DIR}/${run}.fetched" ] && [ ! -e "${STATE_DIR}/${run}.failed" ]; then
                echo "ABORT: the prefetcher exited before ${run} was ready." >&2
                exit 1
            fi
            sleep 15
        done
        waited=$(( SECONDS - waited ))
        fetch_secs="$(cat "${STATE_DIR}/${run}.fetched" 2>/dev/null || echo 0)"
    else
        waited=0
        fetch_start="${SECONDS}"
        fetch_run_pair "${run}" || touch "${STATE_DIR}/${run}.failed"
        fetch_secs=$(( SECONDS - fetch_start ))
    fi
    touch "${STATE_DIR}/${run}.consumed"

    if [ -e "${STATE_DIR}/${run}.failed" ]; then
        echo "  ${run}: retrieval failed, skipping (see ${FETCH_LOG_DIR}/${run}_*.wget.log)" >&2
        FAILED+=("${run}")
        SUMMARY+=("${run}|fetch failed|-|-")
        continue
    fi
    echo "  fetch: ${fetch_secs}s (waited ${waited}s for it)"

    dry="$(snakemake -n "${SNAKEMAKE_ARGS[@]}" "${cram_target}" 2>&1)"
    # The whole point of the seeded manifest and the download is that neither
    # of these rules runs. If one is still in the DAG, continuing would quietly
    # start a multi-hour prefetch instead.
    if printf '%s\n' "${dry}" | grep -qE '^[[:space:]]*(rule|checkpoint) (fetch_hcc1395_fastq|resolve_hcc1395_runs):'; then
        printf '%s\n' "${dry}" >&2
        echo "ABORT: the dry run still schedules fetch_hcc1395_fastq/resolve_hcc1395_runs for ${run}." >&2
        exit 1
    fi

    # One invocation for both bam_target and cram_target, not two: Snakemake
    # resolves the whole DAG up front, so align_hcc1395 and cram_hcc1395 run
    # back to back in the same process with no rebuild in between. Requesting
    # them as two separate snakemake calls (as smoke_hcc1395_local_driver.sh
    # does, deliberately, for a per-stage wall clock) is unsafe here:
    # run_hcc1395_alignment.py deletes its input FASTQ once align_hcc1395
    # succeeds, so a second, later invocation rebuilds the DAG from scratch,
    # sees fetch_hcc1395_fastq's outputs missing, and reschedules it -- through
    # the slow SRA path the whole point of this driver is to avoid -- and then
    # treats align_hcc1395 as stale too, since the re-fetched FASTQ postdates
    # the sorted.bam already sitting there, redoing the multi-hour alignment.
    # Confirmed live on 2026-09-03: exactly that happened on SRR7890858.
    stage_start="${SECONDS}"
    if ! snakemake "${SNAKEMAKE_ARGS[@]}" "${bam_target}" "${cram_target}"; then
        echo "  ${run}: align/cram failed, skipping" >&2
        FAILED+=("${run}"); SUMMARY+=("${run}|${fetch_secs}|align/cram failed|-"); continue
    fi
    align_cram_secs=$(( SECONDS - stage_start ))

    place_secs=0
    if [ "${SKIP_HDD_PLACE:-0}" != "1" ]; then
        stage_start="${SECONDS}"
        mkdir -p "${run_dir}"
        rsync -a --partial "${cram_target}" "${cram_target}.crai" "${run_dir}/"
        if [ "$(stat -c%s "${cram_target}")" != "$(stat -c%s "${run_dir}/${run}.sorted.cram")" ]; then
            echo "ABORT: ${run_dir}/${run}.sorted.cram differs in size from ${cram_target}; keeping the scratch copy." >&2
            exit 1
        fi
        conda run -n flash_fraq samtools quickcheck -v "${run_dir}/${run}.sorted.cram"
        place_secs=$(( SECONDS - stage_start ))
        if [ "${KEEP_SCRATCH:-0}" != "1" ]; then
            rm -f "${cram_target}" "${cram_target}.crai"
        fi
    fi

    SUMMARY+=("${run}|${fetch_secs}|${align_cram_secs}|${place_secs}")
    echo "  done: align+cram ${align_cram_secs}s, place ${place_secs}s"
done

[ -z "${PREFETCH_PID}" ] || wait "${PREFETCH_PID}" 2>/dev/null || true

# --- 7. summary -------------------------------------------------------------

echo
echo "=== ${#PENDING[@]} run(s) on $(hostname) ($(nproc) cores) ==="
printf '  %-12s %10s %14s %10s\n' run fetch align+cram place
for row in "${SUMMARY[@]}"; do
    IFS='|' read -r r f ac p <<< "${row}"
    printf '  %-12s %10s %14s %10s\n' "${r}" "${f}" "${ac}" "${p}"
done
echo "  (fetch seconds are the download's own duration -- for every run after"
echo "   the first it overlapped the previous run's alignment, so it is off the"
echo "   critical path)"

if [ ${#FAILED[@]} -gt 0 ]; then
    echo
    echo "FAILED: ${FAILED[*]}" >&2
fi

echo
echo "Note: '${HDD_DATASET}' is not in local_features.yaml's datasets: list, so these"
echo "CRAMs will not enter the feature pipeline until it is added there."

[ ${#FAILED[@]} -eq 0 ]
