#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G  
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-samtools-%j.out
#SBATCH --job-name=samtools
#SBATCH --time=02:00:00


set -euo pipefail
#Running samtools 

# Container
SAMTOOLS_CONTAINER=oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f

# Directories
bam=$1
BAM_DIR=$2
QC_OUT_DIR=$3
mkdir -p "$QC_OUT_DIR"

echo "# Starting samtools QC using container"
date
echo "# BAM directory: $BAM_DIR"
echo "# Output QC dir: $QC_OUT_DIR"
echo



    #Set sample to be the identifier code (Ex: SRR1555695)
    sample=$(basename "$bam" _Aligned.sortedByCoord.out.bam)

    echo "--------------------------------------------"
    echo "# Processing sample: $sample"
    echo

    # make sure BAM is indexed
    if [[ ! -f "${bam}.bai" ]]; then
    # -f is filename. ! inverts the truth; if ! were absent it would mean "if the bam file index (bai) is present, then ..." with ! it means "if the bam file index bai is absent, then..." 
        apptainer exec "$SAMTOOLS_CONTAINER" samtools index "$bam"
    #This command creates the missing bai index file from the bam file for later commands
    fi

    # flagstat
    apptainer exec "$SAMTOOLS_CONTAINER" samtools flagstat \
        "$bam" > "$QC_OUT_DIR"/"$sample".flagstat.txt
#Gives alignment summary to see if reads aligned in a way that makes sense: 
#- total reads
#- mapped reads
#- properly paired reads
#- singletons
#- secondary or supplementary alignments
#- duplicates
    # stats
    apptainer exec "$SAMTOOLS_CONTAINER" samtools stats \
        "$bam" > "$QC_OUT_DIR"/"$sample".samstats.txt
#Gives text report for use with Multiqc 

    # idxstats (uses the .bai index)
    apptainer exec "$SAMTOOLS_CONTAINER" samtools idxstats \
        "$bam" > "$QC_OUT_DIR"/"$sample".idxstats.txt
#Gives a table of read counts per sequence by chromosome, length, mapped_reads and unmapped_reads to give a distribution 



echo
echo "# samtools QC complete"
apptainer exec "$SAMTOOLS_CONTAINER" samtools --version
date