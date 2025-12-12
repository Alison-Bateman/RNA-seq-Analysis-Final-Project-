#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G  
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-star_align-%j.out
#SBATCH --job-name=star_align_all

set -euo pipefail

# Container
STAR_CONTAINER=oras://community.wave.seqera.io/library/star:2.7.11b--84fcc19fdfab53a4

# Positional arguments:
FASTQ_DIR=$1
INDEX_DIR=$2
GTF_FILE=$3
OUTDIR=$4

#make outdir directory 
mkdir -p "$OUTDIR"

#Initial Logging 
echo "# STAR batch alignment started"
date
echo "# FASTQ directory: $FASTQ_DIR"
echo "# Genome index:   $INDEX_DIR"
echo "# GTF:            $GTF_FILE"
echo "# Output dir:     $OUTDIR"
echo

# Loop through all R1 files
for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"

    # Extract sample ID
    sample_id=$(basename "$R1" _1.fastq.gz)

    echo "--------------------------------------------"
    echo "# Processing sample: $sample_id"
    echo "# R1: $R1"
    echo "# R2: $R2"
    echo

    apptainer exec "$STAR_CONTAINER" STAR \
        --runThreadN "$SLURM_CPUS_PER_TASK" \
        --genomeDir "$INDEX_DIR" \
        --readFilesIn "$R1" "$R2" \
        --readFilesCommand zcat \
        --sjdbGTFfile "$GTF_FILE" \
        --outFileNamePrefix "$OUTDIR"/"${sample_id}_" \
        \
        --outFilterType BySJout \
        --outFilterMultimapNmax 50 \
        --winAnchorMultimapNmax 100 \
        --outSAMmultNmax 10 \
        --outFilterMismatchNoverReadLmax 0.03 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        \
        --quantMode TranscriptomeSAM \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes Standard

# --outFilterType BySJout ==> Filter alignments based on supported splice junctions.
#      STAR keeps only alignments consistent with high-confidence junctions, reducing false-positive
#      splice junctions (helpful in large, repetitive, polyploid genomes like hexaploid Helianthus).

# --outFilterMultimapNmax 50 ==> Maximum number of loci a read is allowed to map to.
#      Reads mapping to >50 places are discarded. This keeps highly repetitive reads from overwhelming
#      the data, but still allows some multi-mapping in a complex genome.

# --winAnchorMultimapNmax 100 ==> Max number of loci for each seed/anchor in the seed-search step.
#      Limits how many candidate locations STAR explores for each seed, controlling runtime and
#      memory in repetitive regions.

# --outSAMmultNmax 10 ==> Maximum number of alignments per read written to the SAM/BAM.
#      Even if STAR finds up to 50 mapping locations, only the 10 best will be saved in the BAM.
#      This keeps file size reasonable while still recording major multi-mapping cases. Especially importnat for
#      a hexploid genome. 

# --outFilterMismatchNoverReadLmax 0.03 ==> Maximum allowed mismatch rate per read (3% of read length).
#      For ~100 bp reads, this allows ~3 mismatches per read. Tight enough to remove low-quality or
#      mismapped reads, but not so strict that natural variation or sequencing errors break mapping. Slightly more 
#      stringent than whats standard. 

# --alignSJoverhangMin 8 ==> Minimum number of bases that must align on each side of a *novel* splice junction.
#      Requires at least 8 bp into each exon, which stabilizes junction calling and reduces incorrect
#      junctions that come from short, low-confidence overhangs.

# --alignSJDBoverhangMin 1 ==> Minimum overhang for *annotated* (GTF-based) splice junctions.
#      Very permissive for known junctions: even 1 bp overhang can support an annotated junction.
#      This favors known GTF junctions while still enforcing stricter rules on novel ones.

# --alignIntronMin 20 ==> Minimum intron length allowed.
#      Filters out tiny “introns” that are likely just small deletions or misalignments. 20 bp is standard.

# --alignIntronMax 1000000 ==> Maximum intron length (1 Mb).
#      Allows very long introns (as seen in plants) but prevents absurdly huge gaps being treated as introns. Standard Value. 

# --alignMatesGapMax 1000000 ==> Maximum allowed genomic distance between paired mates (1 Mb).
#      Prevents mates that align extremely far apart (likely chimeric / wrong) from being treated as a
#      proper pair, but still accommodates very large gene/intergenic regions in big genomes. Standard value. 

# --quantMode TranscriptomeSAM ==> In addition to the genomic BAM, output a transcriptome-aligned BAM
#      (Aligned.toTranscriptome.out.bam). This projects reads onto transcript coordinates using the GTF.
#      This transcriptome BAM is REQUIRED for alignment-based quantification with Salmon for expression estimation / differential expression.

# --outSAMtype BAM SortedByCoordinate ==> Write the main alignment file as a BAM, sorted by genomic coordinate.
#      This is the standard format for downstream tools (samtools, MultiQC, IGV etc.).

# --outSAMattributes Standard ==> Include a standard set of SAM tags (NH, HI, AS, nM, etc.) in the BAM.
#      These tags are useful for QC, multi-mapping information, and for tools that rely on alignment scores.

    echo "# Finished sample: $sample_id"
    date
done

echo
echo "# All samples complete"
apptainer exec "$STAR_CONTAINER" STAR --version
date

