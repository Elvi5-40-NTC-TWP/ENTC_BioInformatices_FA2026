#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ACCESSIONS_FILE="../Accessions.txt"
SRR_DIR="SRR"
PREP_SCRIPT="$SCRIPT_DIR/Elvis-NTC-General-Prep.sh"
POP_PREP_SCRIPT="$SCRIPT_DIR/Elvis-NTC-Population_Analysis_Prep.sh"

mkdir -p "$SRR_DIR" logs

# Read accessions into an array
ACCESSIONS=()
while read -r acc || [ -n "$acc" ]; do
    acc=$(echo "$acc" | tr -d '\r' | xargs)
    [ -z "$acc" ] && continue
    ACCESSIONS+=("$acc")
done < "$ACCESSIONS_FILE"

# Remove any accession that already has both a .sorted.bam and .sorted.bam.bai in bam_out/
FILTERED_ACCESSIONS=()
for acc in "${ACCESSIONS[@]}"; do
    if [ -f "bam_out/${acc}.sorted.bam" ] && [ -f "bam_out/${acc}.sorted.bam.bai" ]; then
        echo "Removing $acc from array (already in bam_out/)"
        continue
    fi
    FILTERED_ACCESSIONS+=("$acc")
done

# Load SRA toolkit if needed
module load sratoolkit

# Prefetch remaining .sra files into the SRR folder first
for acc in "${FILTERED_ACCESSIONS[@]}"; do
    echo "Prefetching $acc into $SRR_DIR..."
    prefetch --output-directory "$SRR_DIR" "$acc"

    sleep 2
done

echo "All prefetches finished. Submitting $PREP_SCRIPT array job..."
sbatch "$PREP_SCRIPT"
echo "$PREP_SCRIPT is done. Time for $POP_PREP_SCRIPT job..."
sbatch "$POP_PREP_SCRIPT"