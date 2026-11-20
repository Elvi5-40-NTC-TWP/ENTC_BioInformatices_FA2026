#!/bin/bash

ACCESSIONS_FILE="Accessions.txt"
SRR_DIR="SRR"
PREP_SCRIPT="Elvis-NTC-General-Prep.sh"

mkdir -p "$SRR_DIR" logs

# Load SRA toolkit if needed
module load sra 2>/dev/null || true

# Prefetch all .sra files into the SRR folder first
while read -r acc || [ -n "$acc" ]; do
    # Skip empty lines or carriage returns
    acc=$(echo "$acc" | tr -d '\r' | xargs)
    [ -z "$acc" ] && continue

    echo "Prefetching $acc into $SRR_DIR..."
    prefetch --output-directory "$SRR_DIR" "$acc"

    sleep 2
done < "$ACCESSIONS_FILE"

echo "All prefetches finished. Submitting $PREP_SCRIPT array job..."
sbatch "$PREP_SCRIPT"