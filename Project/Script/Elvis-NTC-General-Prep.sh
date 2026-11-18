#!/bin/bash

#SBATCH --job-name=ENTC_LABEight_PopGen_Pipeline
#SBATCH --time=03:00:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=5
#SBATCH --partition=batch
#SBATCH --array=1-28%5
#SBATCH --output=logs/sra_job_%A_%a.out
#SBATCH --error=logs/sra_job_%A_%a.errD)