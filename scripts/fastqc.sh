#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-fastqc-%j.out
#SBATCH --job-name=fastqc
#SBATCH --time=01:00:00
#SBATCH --mem=4G

set -euo pipefail

# Load the OSC module for FastQC
module load fastqc/0.12.1

# Positional arguments:
fastq_file=$1
FASTQC_OUT_DIR=$2

# Initial logging
echo "Starting script fastqc.sh"
date
echo "Input FASTQ file:  $fastq_file"
echo "Output dir:         $FASTQC_OUT_DIR"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$FASTQC_OUT_DIR"

# Run FastQC
fastqc \
    --threads 4 \
    --outdir "$FASTQC_OUT_DIR" \
    $fastq_file

#Finished sample logging
echo "Finished fastqc.sh on $fastq_file" 
date

# Final logging
echo "Script complete"