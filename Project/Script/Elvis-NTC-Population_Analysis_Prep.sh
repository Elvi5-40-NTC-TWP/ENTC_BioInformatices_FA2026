#!/bin/bash

#SBATCH --job-name=ENTC_LABEight_PopGen_Pipeline
#SBATCH --output=pipeline_%j.log
#SBATCH --error=pipeline_%j.err
#SBATCH --time=03:00:00
#SBATCH --mem=15G
#SBATCH --cpus-per-task=5
#SBATCH --partition=batch