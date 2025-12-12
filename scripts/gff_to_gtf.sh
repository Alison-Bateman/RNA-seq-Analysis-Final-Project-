#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --output=gff_convert_%j.out
#SBATCH --error=gff_convert_%j.err
#SBATCH --job-name=GFF_to_GTF

set -euo pipefail

# Positional arguments: 
GFF_FILE=$1
GTF_OUT_FILE=$2
OUT_DIR=$3

mkdir -p $OUT_DIR
# initial Logging
echo 
echo "Converting $GFF_FILE to GTF file at $GTF_OUT_FILE"

mkdir -p $OUT_DIR
# Run Command 
apptainer exec oras://community.wave.seqera.io/library/gffread:0.12.7--b08e770b84a4a126 \
gffread $GFF_FILE -T -o $GTF_OUT_FILE

#Final logging
echo "Completed File Conversion"