#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=16
#SBATCH --output=slurm-fastqc-%j.out
#SBATCH --cpus-per-task=8

set -euo pipefail

# Positional arguments: 
RUN=$1
OUT_DIR=$2

#Initial Logging
echo "Starting fastq_download&split.sh"
date
echo "NCBI Run:   $RUN"
echo "Output dir: $OUT_DIR"
fasterq-dump --version
echo "-------------------------------------------------------"

#Make out dir if it doesn't already exist
mkdir -p $OUT_DIR

echo "Starting Prefetch from NBCI for $RUN"
# Download the raw run from NCBI 
 # 1) Download the run using the SRA identifier code (Ex: SRR1555696)
prefetch -O $OUT_DIR $RUN  
echo "Download Complete"
echo "-------------------------------------------------------"

# 2) Convert to FASTQ and split mates
echo "Running fasterq-dump on $RUN"
fasterq-dump $RUN --split-files --threads 16 -O $OUT_DIR

echo "$RUN converted fastq file with split mates"
# yields fastq/SRR1555696_1.fastq and fastq/SRR1555696_2.fastq

echo "Downloaded and split $RUN"
echo "-------------------------------------------------------"

echo "compressing files"
# 3) Compress files 
apptainer exec \
oras://community.wave.seqera.io/library/pbgzip:2016.08.04--c8e0d8eb135a301c \
  pbgzip -n 8 $OUT_DIR/"$RUN"_1.fastq \
  pbgzip -n 8 $OUT_DIR/"$RUN"_2.fastq \

# yields fastq/SRR1555696_1.fastq.gz and fastq/SRR1555696_2.fastq.gz
echo "Files are compressed and you are loved" 
echo "-------------------------------------------------------"
echo "script complete"