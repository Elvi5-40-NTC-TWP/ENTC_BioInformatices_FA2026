#!/bin/bash

#SBATCH --job-name=Elvis-NTC-Population_Analysis_Prep_Pipeline
#SBATCH --output=pipeline_%j.log
#SBATCH --error=pipeline_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=5
#SBATCH --partition=batch

set -euo pipefail

export HOME="${HOME:-/home/entruji2}"
export USER="${USER:-entruji2}"
export PATH="/usr/local/bin:/usr/bin:/bin:${PATH:-}"

module load python
module load sratoolkit
module load java

PROJECT_BASE="/project/stuckert/entruji2/ENTC_BioInformatices_FA2026/Project/White-Space"
BAM_DIR="${PROJECT_BASE}/bam_out"
POPGEN_BASE="${PROJECT_BASE}/popgen_out"
CHR_DIR="${POPGEN_BASE}/chr_beagle"
REF_GENOME="${PROJECT_BASE}/ncbi_dataset/data/GCA_034140825.1/GCA_034140825.1_AGIS1.0_genomic.fna"
BAM_LIST="${POPGEN_BASE}/bam_list.txt"
ACCESSIONS_FILE="${PROJECT_BASE}/../Accessions.txt"
META_TSV="${PROJECT_BASE}/../metadata.tsv"

mkdir -p "${POPGEN_BASE}" "${CHR_DIR}"
cd "${POPGEN_BASE}"

if [ -f "${ACCESSIONS_FILE}" ]; then
    {
        printf "sample_id\tgroup\n"
        while read -r acc || [ -n "$acc" ]; do
            acc=$(echo "$acc" | tr -d '\r' | xargs)
            [ -z "$acc" ] && continue
            [ -f "${BAM_DIR}/${acc}.sorted.bam" ] && [ -f "${BAM_DIR}/${acc}.sorted.bam.bai" ] || continue
            printf "%s\tUNASSIGNED\n" "$acc"
        done < "${ACCESSIONS_FILE}"
    } > "${META_TSV}"
fi

echo "=== Step 0.5: Re-indexing Reference FASTA if Needed ==="
if [ ! -f "${REF_GENOME}.fai" ] || [ "${REF_GENOME}" -nt "${REF_GENOME}.fai" ]; then
    samtools faidx "${REF_GENOME}"
fi

echo "=== Step 1: Generating BAM list ==="
find "${BAM_DIR}" -maxdepth 1 -name "*.sorted.bam" | sort > "${BAM_LIST}"

echo "=== Step 2: Running ANGSD Per-Chromosome (-r) for Beagle GLs ==="
COMBINED_BEAGLE="${POPGEN_BASE}/angsd_output.beagle.gz"
rm -f "${COMBINED_BEAGLE}" "${CHR_DIR}"/*.beagle.gz

FIRST_CHR=1
while IFS=$'\t' read -r chrom _rest; do
    [ -z "$chrom" ] && continue
    chr_prefix="${CHR_DIR}/angsd_${chrom}"

    angsd -GL 1 \
        -r "${chrom}" \
        -out "${chr_prefix}" \
        -nThreads 5 \
        -ref "${REF_GENOME}" \
        -bam "${BAM_LIST}" \
        -minMapQ 30 -minQ 20 \
        -doMajorMinor 1 \
        -SNP_pval 1e-6 -doMaf 1 \
        -minMaf 0.05 \
        -skipTriallelic 1 \
        -doGlf 2

    chr_beagle="${chr_prefix}.beagle.gz"
    if [ -f "${chr_beagle}" ]; then
        if [ "${FIRST_CHR}" -eq 1 ]; then
            zcat "${chr_beagle}" | gzip -c > "${COMBINED_BEAGLE}"
            FIRST_CHR=0
        else
            zcat "${chr_beagle}" | tail -n +2 | gzip -c >> "${COMBINED_BEAGLE}"
        fi
    fi
done < "${REF_GENOME}.fai"

echo "=== Step 3: Running PCAngsd ==="
pcangsd -b "${COMBINED_BEAGLE}" \
    --admix \
    --snp-weights \
    --sites-save \
    --threads 5 \
    --out "${POPGEN_BASE}/pcangsd_output"

echo "Pipeline execution completed successfully!"