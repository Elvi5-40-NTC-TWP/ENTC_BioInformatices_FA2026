#!/bin/bash
#SBATCH --job-name=Elvis-NTC-General-Prep_Pipeline
#SBATCH --time=03:00:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=5
#SBATCH --partition=batch
#SBATCH --array=1-28%2
#SBATCH --output=logs/sra_job_%A_%a.out
#SBATCH --error=logs/sra_job_%A_%a.err

set -euo pipefail

# Ensure .bashrc mappings (Labs and subfolders) and required tools are loaded
if [ -f "$HOME/.bashrc" ]; then
    set +u
    source "$HOME/.bashrc"
    set -u
fi

module load sratoolkit
module load java
module load python

mkdir -p logs fastq_raw fastq_trimmed fastqc_out bam_out bwa_idx tmp_work

TMP_DIR="tmp_work/fasterq_${SLURM_JOB_ID:-$$}_${SLURM_ARRAY_TASK_ID:-0}"
mkdir -p "$TMP_DIR"
trap 'rm -rf "$TMP_DIR"' EXIT

ACCESSIONS_FILE="../Accessions.txt"
SRR_DIR="SRR"
REF_FASTA="ncbi_dataset/data/GCA_034140825.1/GCA_034140825.1_AGIS1.0_genomic.fna"
BWA_EXE="/project/stuckert/entruji2/Labs/Lab_Six/bwa-mem2-2.2.1_x64-linux/bwa-mem2"
BWA_PREFIX="bwa_idx/agis"

# Grab the accession for this array task ID (or optional $1 argument if tested manually)
if [ -n "${1:-}" ]; then
    SAMPLE="$1"
else
    SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCESSIONS_FILE" | tr -d '\r' | xargs)
fi

if [ -z "$SAMPLE" ]; then
    echo "Error: No sample accession found for task ${SLURM_ARRAY_TASK_ID:-manual}" >&2
    exit 1
fi

echo "=== Processing $SAMPLE (Task ${SLURM_ARRAY_TASK_ID:-manual}) ==="

# Locate the prefetched .sra file/directory
SRA_INPUT="$SRR_DIR/${SAMPLE}/${SAMPLE}.sra"
if [ ! -f "$SRA_INPUT" ]; then
    SRA_INPUT="$SRR_DIR/${SAMPLE}.sra"
fi

if [ ! -f "$SRA_INPUT" ]; then
    echo "Error: SRA file not found for $SAMPLE in $SRR_DIR" >&2
    exit 1
fi

# 1. fasterq-dump into uncompressed fastqs, immediately gzip, and remove raw .sra
fasterq-dump --split-files -e 5 -t "$TMP_DIR" --outdir fastq_raw "$SRA_INPUT"
gzip -f "fastq_raw/${SAMPLE}_1.fastq" "fastq_raw/${SAMPLE}_2.fastq"
rm -rf "$SRR_DIR/${SAMPLE}" "$SRR_DIR/${SAMPLE}.sra"

# 2. FastQC on raw reads (--outdir rather than -o)
fastqc -t 5 --outdir fastqc_out "fastq_raw/${SAMPLE}_1.fastq.gz" "fastq_raw/${SAMPLE}_2.fastq.gz"

# 3. fastp trimming & filtering, then remove raw fastq.gz to save space
fastp -w 5 \
  -i "fastq_raw/${SAMPLE}_1.fastq.gz" \
  -I "fastq_raw/${SAMPLE}_2.fastq.gz" \
  -o "fastq_trimmed/${SAMPLE}_1.trim.fastq.gz" \
  -O "fastq_trimmed/${SAMPLE}_2.trim.fastq.gz" \
  --html "fastqc_out/${SAMPLE}_fastp.html" \
  --json "fastqc_out/${SAMPLE}_fastp.json"

rm -f "fastq_raw/${SAMPLE}_1.fastq.gz" "fastq_raw/${SAMPLE}_2.fastq.gz"

# 4. Build bwa-mem2 index once with short prefix (protected with mkdir lock for array safety)
if [ ! -f "${BWA_PREFIX}.bwt.2bit.64" ]; then
    if mkdir "${BWA_PREFIX}.lock" 2>/dev/null; then
        "$BWA_EXE" index -p "$BWA_PREFIX" "$REF_FASTA"
        rmdir "${BWA_PREFIX}.lock"
    else
        while [ -d "${BWA_PREFIX}.lock" ]; do
            sleep 10
        done
    fi
fi

# 5. bwa-mem2 alignment -> sorted BAM -> indexed BAM, then clean up trimmed FASTQs
"$BWA_EXE" mem -t 5 "$BWA_PREFIX" \
  "fastq_trimmed/${SAMPLE}_1.trim.fastq.gz" \
  "fastq_trimmed/${SAMPLE}_2.trim.fastq.gz" \
  | samtools sort -@ 5 -m 1500M -o "bam_out/${SAMPLE}.sorted.bam" -

samtools index -@ 5 "bam_out/${SAMPLE}.sorted.bam"

rm -f "fastq_trimmed/${SAMPLE}_1.trim.fastq.gz" "fastq_trimmed/${SAMPLE}_2.trim.fastq.gz"

echo "=== Finished $SAMPLE ==="