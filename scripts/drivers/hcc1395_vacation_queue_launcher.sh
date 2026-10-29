#!/usr/bin/env bash
# Chain the long unattended HCC1395 queue behind the in-flight matched-pair run.
#
# Waits for hcc1395_matched_pair_local_driver.sh to exit (it holds the same
# flock the ENA driver takes, so starting early would just block or collide),
# confirms the ladder-unblocking pair actually produced two CRAMs, then works
# down manifests/hcc1395/hcc1395_vacation_queue.txt.
#
# The queue is deliberately LONGER than the time available: it is ordered by
# value (cross-site T/N pairs, then SPP pooled mixtures, then extra lanes) and
# simply gets as far as it gets. Runs are capped at 220 Gbp because
# src/flash_fraq/alignment.py writes an unsorted BAM to disk before sorting to
# a second file, so peak scratch is FASTQ + unsorted + sorted ~= 1.81 GB/Gbp;
# larger runs cannot fit the 660 GB NVMe and would fail hours in, on a disk
# that the next run then also fails on.
#
# PREFETCH_AHEAD=1 (not the default 2): with runs up to 220 Gbp, two staged
# FASTQ pairs plus an active run's 398 GB peak overruns the disk. The driver's
# MIN_FREE_GB floor gates the prefetcher dynamically as a second line of
# defence.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

QUEUE_FILE="manifests/hcc1395/hcc1395_vacation_queue.txt"
CRAM_ROOT="/mnt/8tb_hdd/flash_fraq_crams/HCC1395_SEQC2/CRAM"
WAIT_PID="${WAIT_PID:-}"
LOG="logs/hcc1395_vacation_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee -a "${LOG}") 2>&1
echo "launcher started $(date -Is); log ${LOG}"

if [ -n "${WAIT_PID}" ]; then
    echo "waiting for matched-pair driver (pid ${WAIT_PID}) to finish..."
    while kill -0 "${WAIT_PID}" 2>/dev/null; do sleep 60; done
    echo "matched-pair driver exited $(date -Is)"
fi

# Do not start a 9-day queue if the pair it is supposed to build on failed.
missing=0
for run in SRR7890858 SRR7890863; do
    if [ -s "${CRAM_ROOT}/${run}/${run}.sorted.cram" ]; then
        echo "  pair OK: ${run}"
    else
        echo "  pair MISSING: ${run} has no CRAM" >&2
        missing=1
    fi
done
if [ "${missing}" -ne 0 ]; then
    echo "ABORT: the matched pair did not complete; not starting the queue." >&2
    exit 1
fi

mapfile -t RUNS < <(grep -vE '^\s*(#|$)' "${QUEUE_FILE}")
echo "queue: ${#RUNS[@]} runs"

# Three passes, not one. A run whose fetch exhausts its wget retries is marked
# <run>.failed and skipped for the rest of that pass -- permanently, as far as
# one invocation is concerned. MedUni Graz IT announced short network outages
# on 2026-09-07/08/10 between 20:00 and 24:00, so a fetch can die for reasons
# that have nothing to do with the run being bad.
#
# Re-running the whole queue is cheap and safe: already_done() in the driver
# skips any run whose CRAM and .crai are already on the HDD, so a later pass
# only re-attempts what actually failed. The hour between passes lets a
# transient outage clear rather than burning the retry on the same dead link.
for pass in 1 2 3; do
    echo
    echo "######## pass ${pass} of 3 -- $(date -Is) ########"
    if PREFETCH_AHEAD=1 ./scripts/drivers/hcc1395_ena_local_driver.sh "${RUNS[@]}"; then
        echo "pass ${pass}: every queued run is on the HDD; nothing left to retry"
        break
    fi
    echo "pass ${pass}: finished with failures (see FAILED above)" >&2
    # A plain `[ ... ] && { ... }` here would abort the script under `set -e` on
    # pass 3, where the test is false and it is the loop body's last command --
    # a three-pass failure would then end silently, which is the one outcome
    # that most needs to be legible on return.
    if [ "${pass}" -lt 3 ]; then
        echo "sleeping 1h before the next pass"
        sleep 3600
    fi
done
echo
echo "launcher finished $(date -Is)"
