#!/bin/bash
#SBATCH --job-name=Elvis-NTC-General-Prep_Pipeline
#SBATCH --time=03:00:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=5
#SBATCH --partition=batch
#SBATCH --array=1-28%5
#SBATCH --output=logs/sra_job_%A_%a.out
#SBATCH --error=logs/sra_job_%A_%a.err

mkdir -p logs fastq_raw fastq_trimmed fastqc_out bam_out bwa_index

ACCESSIONS_FILE="Accessions.txt"
SRR_DIR="SRR"
REF_FASTA="ncbi_dataset/data/GCA_034140825.1/GCA_034140825.1_AGIS1.0_genomic.fna"
BWA_PREFIX="bwa_index/AGIS1.0_bwa"

# Grab the accession for this array task ID (or optional $1 argument if tested manually)
if [ -n "$1" ]; then
    SAMPLE="$1"
else
    SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCESSIONS_FILE" | tr -d '\r' | xargs)
fi

if [ -z "$SAMPLE" ]; then
    echo "Error: No sample accession found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "=== Processing $SAMPLE (Task ${SLURM_ARRAY_TASK_ID:-manual}) ==="

# Locate the prefetched .sra file/directory
SRA_INPUT="$SRR_DIR/${SAMPLE}/${SAMPLE}.sra"
if [ ! -f "$SRA_INPUT" ]; then
    SRA_INPUT="$SRR_DIR/${SAMPLE}.sra"
fi

# 1. fasterq-dump into uncompressed fastqs, immediately gzip, and remove raw .sra
fasterq-dump --split-files -e 5 -t /tmp --outdir fastq_raw "$SRA_INPUT"
gzip -f "fastq_raw/${SAMPLE}_1.fastq" "fastq_raw/${SAMPLE}_2.fastq"
rm -rf "$SRR_DIR/${SAMPLE}" "$SRR_DIR/${SAMPLE}.sra"

# 2. FastQC on raw reads
fastqc -t 5 -o fastqc_out "fastq_raw/${SAMPLE}_1.fastq.gz" "fastq_raw/${SAMPLE}_2.fastq.gz"

# 3. fastp trimming & filtering, then remove raw fastq.gz to save space
fastp -w 5 \
  -i "fastq_raw/${SAMPLE}_1.fastq.gz" \
  -I "fastq_raw/${SAMPLE}_2.fastq.gz" \
  -o "fastq_trimmed/${SAMPLE}_1.trim.fastq.gz" \
  -O "fastq_trimmed/${SAMPLE}_2.trim.fastq.gz" \
  --html "fastqc_out/${SAMPLE}_fastp.html" \
  --json "fastqc_out/${SAMPLE}_fastp.json"

rm -f "fastq_raw/${SAMPLE}_1.fastq.gz" "fastq_raw/${SAMPLE}_2.fastq.gz"

# 4. Build bwa-mem2 index once (protected with mkdir lock for array safety)
if [ ! -f "${BWA_PREFIX}.bwt.2bit.64" ]; then
    if mkdir "${BWA_PREFIX}.lock" 2>/dev/null; then
        bwa-mem2 index -p "$BWA_PREFIX" "$REF_FASTA"
        rmdir "${BWA_PREFIX}.lock"
    else
        while [ -d "${BWA_PREFIX}.lock" ]; do
            sleep 10
        done
    fi
fi

# 5. bwa-mem2 alignment -> sorted BAM -> indexed BAM, then clean up trimmed FASTQs
bwa-mem2 mem -t 5 "$BWA_PREFIX" \
  "fastq_trimmed/${SAMPLE}_1.trim.fastq.gz" \
  "fastq_trimmed/${SAMPLE}_2.trim.fastq.gz" \
  | samtools sort -@ 5 -m 1500M -o "bam_out/${SAMPLE}.sorted.bam" -

samtools index -@ 5 "bam_out/${SAMPLE}.sorted.bam"

rm -f "fastq_trimmed/${SAMPLE}_1.trim.fastq.gz" "fastq_trimmed/${SAMPLE}_2.trim.fastq.gz"

echo "=== Finished $SAMPLE ==="