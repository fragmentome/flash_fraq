#!/bin/bash
#SBATCH --job-name=fetch-hg19-reference
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#
# Downloads the UCSC hg19 reference FASTA (chr-prefixed contigs, matching
# the "hg19" build name used throughout docs/phase2-preprocessing-plan.md)
# and indexes it for samtools/bwa so MAT-19/MAT-20 alignment can use it.
# Submitted via scripts/drivers/fetch_hg19_reference_driver.sh.
#
# NOTE: the exact reference FASTA variant (decoy set, chr-prefix) behind the
# existing Snyder/Cristiano hg19 BAMs is still unconfirmed (phase2 §0) --
# after this runs, compare this FASTA's indexed contigs against those BAMs
# via FastqAligner.verify_matches_reference_header before trusting it as the
# alignment target for the full sample sets.

set -euo pipefail

: "${REF_DIR:?}"
: "${REF_FASTA:?}"

mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

wget -q -O hg19.fa.gz "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz"
gunzip -f hg19.fa.gz
if [[ "${REF_DIR}/hg19.fa" != "${REF_FASTA}" ]]; then
    mv hg19.fa "${REF_FASTA}"
fi

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate flash_fraq

samtools faidx "${REF_FASTA}"
bwa index "${REF_FASTA}"

missing=()
for f in "${REF_FASTA}" "${REF_FASTA}.fai" "${REF_FASTA}.amb" "${REF_FASTA}.ann" "${REF_FASTA}.bwt" "${REF_FASTA}.pac" "${REF_FASTA}.sa"; do
    if [[ ! -s "${f}" ]]; then
        missing+=("${f}")
    fi
done

if [[ ${#missing[@]} -gt 0 ]]; then
    echo "ERROR: missing or empty output file(s): ${missing[*]}" >&2
    exit 1
fi
