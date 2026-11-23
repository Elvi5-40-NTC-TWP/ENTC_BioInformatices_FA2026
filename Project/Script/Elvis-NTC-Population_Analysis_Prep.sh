#!/bin/bash

#SBATCH --job-name=Elvis-NTC-Population_Analysis_Prep_Pipeline
#SBATCH --output=pipeline_%j.log
#SBATCH --error=pipeline_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=5
#SBATCH --partition=batch

set -euo pipefail

# Ensure .bashrc mappings (Labs and subfolders) and required tools are loaded
if [ -f "$HOME/.bashrc" ]; then
    set +u
    source "$HOME/.bashrc"
    set -u
fi

module load python
module load sratoolkit
module load java

# 0. Define Directories and Variables
PROJECT_BASE="/project/stuckert/entruji2/ENTC_BioInformatices_FA2026/Project/White-Space"
BAM_DIR="${PROJECT_BASE}/bam_out"
POPGEN_BASE="${PROJECT_BASE}/popgen_out"
REF_GENOME="${PROJECT_BASE}/ncbi_dataset/data/GCA_034140825.1/GCA_034140825.1_AGIS1.0_genomic.fna"
BAM_LIST="${POPGEN_BASE}/bam_list.txt"
ACCESSIONS_FILE="${PROJECT_BASE}/../Accessions.txt"
META_TSV="${PROJECT_BASE}/../metadata.tsv"

mkdir -p "${POPGEN_BASE}"
cd "${POPGEN_BASE}"

# Create starter metadata.tsv from Accessions.txt for accessions with an existing .sorted.bam and .bai
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

# 1. ANGSD Genotype Calling & Likelihoods
echo "=== Step 2: Running ANGSD Genotype Calling and GL Generation ==="
angsd -GL 1 \
    -out "${POPGEN_BASE}/angsd_output" \
    -nThreads 5 \
    -ref "${REF_GENOME}" \
    -bam "${BAM_LIST}" \
    -minMapQ 30 -minQ 20 \
    -doMajorMinor 1 \
    -SNP_pval 1e-6 -doMaf 1 \
    -minMaf 0.05 \
    -skipTriallelic 1 \
    -doGlf 2 \

# 2. Population Structure Analysis (PCAngsd)
echo "=== Step 3: Running PCAngsd ==="
# Note: ANGSD outputs the beagle file as angsd_output.beagle.gz
pcangsd -b "${POPGEN_BASE}/angsd_output.beagle.gz" \
    --admix \
    --snp-weights \
    --sites-save \
    --threads 5 \
    --out "${POPGEN_BASE}/pcangsd_output"

echo "Pipeline execution completed successfully!"