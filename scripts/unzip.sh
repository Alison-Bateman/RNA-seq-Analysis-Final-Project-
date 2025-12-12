#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --job-name=zip
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=8G


set -euo pipefail

for fastq in fastq/*.fastq; do
  apptainer exec \
    oras://community.wave.seqera.io/library/pbgzip:2016.08.04--c8e0d8eb135a301c \
    pbgzip -n 8 "$fastq"
done

echo "zipped"
