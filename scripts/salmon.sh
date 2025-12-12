#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --job-name=salmon_quant_jung
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=slurm-salmon-quant-%j.out

set -euo pipefail

SALMON_CONTAINER=oras://community.wave.seqera.io/library/salmon:1.10.3--726401738a281398

bam=$1
TRANSCRIPTOME_DIR=$2
GTF=$3
SALMON_OUT_DIR=$4


# Make out directory 
mkdir -p "$SALMON_OUT_DIR"

echo

  bamfile=$(basename "$bam")

  # Extract SRR ID (everything before first underscore)
  sample="${bamfile%%_*}"

  outdir="$SALMON_OUT_DIR/$sample"

  echo "==============================================="
  echo "# Sample:     $sample"
  echo "# Input BAM:  $bam"
  echo "# Output dir: $outdir"
  echo "==============================================="

  mkdir -p "$outdir"

  apptainer exec "$SALMON_CONTAINER" salmon quant \
    -t "$TRANSCRIPTOME_DIR" \
    -l A \
    -a "$bam" \
    -p 8 \
    --gcBias \
    --seqBias \
    -g "$GTF" \
    -o "$outdir"
    
# -t ==> input transcriptome 
# -i "$INDEX_DIR" ==> set path to index 
# -l A ==> auto detect library type (strandedness)
# -a "$bam" ==> Alignment-based mode using STAR BAMs
# -p 8 ==> threads per task
# --validateMappings ==> validation for alignments
# --gcBias ==> correct gc bias in reads
# --seqBias ==> corrects sequence-specific bias
# -o "$outdir"



  echo "# Done: $sample"
  echo

