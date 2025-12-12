#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-star_index-%j.out
#SBATCH --mem=160G
#SBATCH --time=24:00:00
set -euo pipefail

# Container
STAR_CONTAINER=oras://community.wave.seqera.io/library/star:2.7.11b--84fcc19fdfab53a4

# Positional arguments:
FASTA_FILE=$1
GTF_FILE=$2
OUTDIR=$3

# Initial logging
echo "# Starting script star_index.sh"
date
echo "# Input assembly FASTA file:      $FASTA_FILE"
echo "# Input annotation GTF file:      $GTF_FILE"
echo "# Output dir:                     $OUTDIR"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p $OUTDIR

# Run STAR
apptainer exec "$STAR_CONTAINER" STAR \
    --runMode genomeGenerate \
    --genomeFastaFiles "$FASTA_FILE" \
    --sjdbGTFfile "$GTF_FILE" \
    --sjdbOverhang 99 \
    --genomeSAindexNbases 14 \
    --genomeDir "$OUTDIR" \
    --runThreadN 16

# Explanation of  options:
#   --runMode genomeGenerate
#       Build a genome index instead of aligning reads.
#   --runThreadN 16
#       Match Slurm's --cpus-per-task (STAR uses 16 threads).
#   --sjdbOverhang
#       Set to max read length - 1. For ~100 bp reads, 99.
#   --genomeSAindexNbases
#       Controls suffix array depth (speed vs RAM).
#       Default ≈ min(14, log2(genome_length)/2 - 1).
#       For ~10.5 Gb genome, that gives ~15.5 → capped to 14,
#       but 14 can spike RAM on huge plant genomes, so we use 13
#       to reduce memory pressure with minimal performance cost.
# Final logging
echo
echo "# Used STAR version:"
apptainer exec "$STAR_CONTAINER" STAR --version
echo "# Successfully finished script star_index.sh"
date
