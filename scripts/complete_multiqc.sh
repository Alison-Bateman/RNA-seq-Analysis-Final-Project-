#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-multiqc-%j.out
#SBATCH --job-name=multiqc

set -euo pipefail


# Constants
MULTIQC_CONTAINER=oras://community.wave.seqera.io/library/multiqc:1.31--e111a2e953475c51

#postional Parameters: 
RAW_FASTQC_DIR=$1
TRIMMED_FASTQC_FASTP_DIR=$2
STAR_RESULTS_DIR=$3
SAMTOOLS_RESULTS_DIR=$4
SALMON_RESULTS_DIR=$5
OUTDIR=$6

#%date for file name 
run_date=$(date +%m%d)

# Initial logging
echo "# Starting script multiqc.sh"
date
echo "# Inputs:                         raw fastqc, trimmed fastqc, STAR results, samtools results"
echo "# Output dir:                     $OUTDIR"
echo

# Create the output dir
mkdir -p "$OUTDIR"

# Run MultiQC
apptainer exec "$MULTIQC_CONTAINER" multiqc \
    --outdir "$OUTDIR" \
    "$RAW_FASTQC_DIR" \
    "$TRIMMED_FASTQC_FASTP_DIR" \
    "$STAR_RESULTS_DIR" \
    "$SAMTOOLS_RESULTS_DIR" \
    "$SALMON_RESULTS_DIR" \
    --filename "full_multiqc_${run_date}.html"

# Final logging
echo "# Finished multiqc script"
date
