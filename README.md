# Project Goal

I'll be going through the RNA seq analysis workflow using FASTQ raw reads from *Helianthus tuberosus* (sunchoke) from a study in 2014, and a sunchoke reference genome from last year(Wang et al., 2024). Because I need to do RNA sequencing in the future with sunchoke, having scripts prepared for me to feed files into for analysis will help streamline my future research.

Wang S, Wang A, Chen R, Xu D, Wang H, Jiang F, Liu H, Qian W, Fan W. Haplotype-resolved chromosome-level genome of hexaploid Jerusalem artichoke provides insights into its origin, evolution, and inulin metabolism. Plant Commun. 2024 Mar 11;5(3):100767. doi: 10.1016/j.xplc.2023.100767. Epub 2023 Nov 17. PMID: 37974403; PMCID: PMC10943552.

Jung WY, Lee SS, Kim CW, Kim H-S, Min SR, Moon JS, et al. (2014) RNA-Seq Analysis and De Novo Transcriptome Assembly of Jerusalem Artichoke (Helianthus tuberosus Linne). PLoS ONE 9(11): e111982. https://doi.org/10.1371/journal.pone.0111982

# **Preliminary steps:**

Set directory:

```         
cd /fs/ess/PAS2880/users/bateman139/final_project
```

Make git ignore:

```         
touch .gitignore 
```

the .gitignore file will house the fastq files and ref files.

Make Starting Directories:

```{bash}
mkdir -p docs fastq ref/figshare;gtf results scripts
```

Copy scripts into scripts directory:

``` bash
cp /fs/ess/PAS2880/users/bateman139/project/scripts/* /fs/ess/PAS2880/users/bateman139/final_project/scripts/
```

Change permissions to execute

```         
chmod u+x scripts/*
```

Download toolkits for downloading and splitting fastq files including NBCI command line tool, Datasets, and an SRA Toolkit using the following command:

```{bash}
#Download Datasets command-line tools
echo "Downloading NCBI Datasets Command-line tools"
#Installed using guide at https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/
#Download datasets: 
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
#Make it executable: 
chmod +x datasets 
echo "Command-line tool download complete"

# Download a prebuilt SRA Toolkit
echo "Starting Toolkit Download"
curl -L -o sratoolkit.tar.gz \
  https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.2.1/sratoolkit.3.2.1-ubuntu64.tar.gz
# Unpack the archive
tar -xzf sratoolkit.tar.gz
tooldir=$(find . -maxdepth 1 -type d -name "sratoolkit.*")
export PATH="$PWD/$tooldir/bin:$PATH"

#To prevent having to reexport everytime use: 
export PATH="/fs/scratch/PAS2880/sratoolkit.3.2.1-ubuntu64/bin:$PATH"
echo "Successfully downloaded toolkit"
```

Datasets will be used to retrieve the fastq files from NBCI. This command can also be used to retrieve other types of data, like the reference genome, however, the reference genome I will be using is only available with annotations from Figshare.

The SRA toolkit contains the function, fasterq-dump, which can be used to split interleaved SRA data into split mates, creating R1_fastq and R2_fastq files.

# Downloading Files

## FASTQ files:

My FASTQ file set is from Jung et al., 2014 (Bioproject PRJNA258432). I will be downloading the RNA seq files from 5 different samples, from different plant tissue of Purple Jerusalem Artichoke:

| Sample/Run ID | Tissue sourced from           |
|---------------|-------------------------------|
| SRR1555737    | Tuber2 sample (mature stage)  |
| SRR1555736    | Tuber1 sample (initial stage) |
| SRR1555697    | Stem sample                   |
| SRR1555695    | Root sample                   |
| SRR1555696    | Leaf sample                   |

The files are interleaved, so they download as a single file containing both paired ends. We will be downloading the splitting the files in one script, the `fastq_download&split.sh` script.

This script first retrieves the desired run from NBCI based on the SRA ID using the `prefetch` command from Datasets to the desired directory. Once downloaded, the file is converted to fastq format with split mates using `fasterq-dump --split-file` from the SRA toolkit. Once the files are converted and split, they're compressed with `gzip fastq/*.fastq`, to compress the new input files.

Because there was more than one run, a list of SRA runs from NBCI was input into runs.txt, and a a loop was created for the values in that file:

```         
echo "SRR1555736
SRR1555737
SRR1555695
SRR1555696
SRR1555697" > docs/runs.txt \
```

```{bash}
RUNS_TXT=$(cat docs/runs.txt)
for RUN in $RUNS_TXT; do
sbatch scripts/fastq_download_and_split.sh "$RUN" fastq
done
```

To check on the batch job progress use:

```{bash}
squeue -u $USER -l
```

## Reference Genome:

Because the reference genome and annotation data were not avaliable on NCBI, the FASTA (.fna) and GFF files were downloaded from figshare (https://figshare.com/articles/dataset/Annotated_reference_genome_of_Helianthus_tuberosus/22491205) using the following commands:

Download Script to use with FASTA and GTF files using the `ref_download.sh` script, which contains the following wget command.

```{bash}
LINK=$1
OUTDIR=$2
FILE_NAME=$3
mkdir -p $OUTDIR/$FILE_NAME
wget \
  --content-disposition \
  --trust-server-names \
  --user-agent="Mozilla/5.0" \
  $LINK \
  -O $OUTDIR/$FILE_NAME
ls -lh $OUTDIR/$FILE_NAME
```

### For the Fasta file, run:

```{bash}
sbatch scripts/ref_download.sh https://figshare.com/ndownloader/files/42806323 ref/figshare HelTub_1.0.fn.gz
```

### For GFF file, run:

```{bash}
sbatch scripts/ref_download.sh https://figshare.com/ndownloader/files/42771826 ref/figshare HelTub_1.0.gff.gz
```

#### Unzip with:

```{bash}
gunzip -c ref/figshare/HelTub_1.0.fn.gz > ref/figshare/HelTub_1.0.fn
gunzip -c ref/figshare/HelTub_1.0.gff.gz > ref/figshare/HelTub_1.0.gff
```

### GFF to GTF File Conversion:

In order to use the annotation data for alignment, it needs to be in a GTF file format. To convert the GFF to GTF file, use the following command within the `gff_to_gtf.sh` script:

```{bash}
sbatch scripts/gff_to_gtf.sh ref/figshare/HelTub_1.0.gff ref/gtf/HelTub_1.0.gtf ref/gtf
```

Directly check the file or validate with:

```{bash}
head ref/gtf/HelTub_1.0.gtf
```

Correct Output:

```{bash}
Htub.Chr01.H5 AUGUSTUS transcript 11603 16467 0.04 + . transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144" Htub.Chr01.H5 AUGUSTUS exon 11603 12143 0.2 + . transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"; Htub.Chr01.H5 AUGUSTUS exon 13567 13766 0.23 + . transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144";
```

If you need to zip the files, use the following command (I used pbgzip and gzip because they're large and take awhile):

``` bash
apptainer exec \
oras://community.wave.seqera.io/library/pbgzip:2016.08.04--c8e0d8eb135a301c \
  pbgzip -n 8 ref/figshare/*
  gzip -1 ref/gtf/*

echo "zipped"
```

# **Trimming**

Preform QC on fastq files (R1 and R2) using FASTQC and MultiQC scripts to see the quality of the files and determine how they should be trimmed.

## Initial QC on FASTQ files

1.  Loop all files located in the fastq dir with `fastqc.sh`

    ```{bash}
    for fastq_file in fastq/*fastq.gz; do
      sbatch scripts/fastqc.sh "$fastq_file" results/1_raw_fastqc
    done

    ```

2.  Next, run Multiqc on the fastqc file dir and set the output for the Mutliqc report.

```{bash}
MULTIQC_CONTAINER=oras://community.wave.seqera.io/library/multiqc:1.31--e111a2e953475c51
apptainer exec "$MULTIQC_CONTAINER" multiqc \
results/1_raw_fastqc \
--outdir results/1.5_raw_multiqc
```

## **Analysis of the MultiQC**

The areas that had Fails or Warnings were:

-   Per Tile Sequence Quality (all read 1s are fails and all read 2s are warnings)

-   Per Base Sequence Content (Either fail or warning)

-   Per Sequence GC Content has a warning for SRR1555737 (1 and 2)

-   Per Base N Content (Warnings for all read 1s)

-   Sequence Duplication Levels (Warn or fail)

-   Over represented seqeunces (present in 737 1 and 2, and 736 1 )

The good things about our reads is that:

-   The average base and sequence length are optimal for alignment.

-   They don't need heavy adapter trimming.

Based on this we need to:

1.  Trim low-quality 3' tails (standard).

2.  Keep reads over 50 bp after trimming for alignment.

3.  Filter reads with high N content.

## Trimming The FASTQ Files

1.  To do this I used fastp. I chose fastp as most of the information, tutorials and guides on Multiqc analysis and trimming used fastp, and I found it easier to pick up.

    Use the loop to run the script `fastp.sh` to trim the reads located in the fastq dir:

```{bash}
for R1 in fastq/*_1.fastq.gz; do
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    sbatch scripts/fastp.sh "$R1" "$R2" results/2_trimming 
done
```

## Validate Trimming (QC)

To run the fastqc script on all the trimmed reads, use the loop:

```{bash}
for fastq_file in results/2_trimming/*fastq.gz; do
sbatch scripts/fastqc.sh $fastq_file results/2.5_trimmed_qc/fastqc
done
```

To run Multiqc, I used the command:

```{bash}
MULTIQC_CONTAINER=oras://community.wave.seqera.io/library/multiqc:1.31--e111a2e953475c51
apptainer exec "$MULTIQC_CONTAINER" multiqc \
results/2.5_trimmed_qc/fastqc \
--outdir results/2.5_trimmed_qc/multiqc
```

On the MultiQC report, the reads looked cleaning and contained less duplications and a lower N content.

On the report, there still appears to be a large N content jump at the 35 read postion, but because we filtered all N reads over 5, this must be a small segment under 5, that STAR should be able to handle.

# **Alignment:**

## STAR Index

In order to run the star alignment script, we first need to create an index using STAR's `genomeGenerate` option. This function parses the genome FASTA and integrates the annotation file to include exon-intron structure, building a suffix-array-based index optimized for alignment.

Prior to running the index script, star-index.sh, ensure that the files are unzipped and validated.

To create the index, we need to input a .fna and gtf file. Our FASTA file ends with .fn, not .fna, which will cause issues. Rename the file to have the .fna file ending before creating the index:

```         
mv ref/figshare/HelTub_1.0.fn ref/figshare/HelTub_1.0.fna
```

Next, we can run the `star_ index.sh` script with the following command, inputting the .fna file, the .gtf file, and the output dir.

**WARNING: THIS BATCH JOB TAKES HOURS TO COMPLETE**

```{bash}
sbatch scripts/star_index.sh ref/figshare/HelTub_1.0.fna ref/gtf/HelTub_1.0.gtf ref/index
```

## STAR Alignment

Run the `star_align.sh` script with the following command, inputting the fastq dir, the index dir, the gtf file dir, and the results out dir.

```{bash}
sbatch scripts/star_align.sh fastq ref/index ref/gtf/HelTub_1.0.gtf results/3_star
```

# SamTools

Samtools lets us inspect, summarize and manage data in BAM files. The script I created runs four different commands under Samtools; index, flagstat, stats, and idxstats. Index creates BAI files from existing BAM files created by STAR alignment, to be used for the idxstats command. flagstats gives an alignment summary to see if reads aligned in a way that makes sense, including the total reads, mapped reads, properly paired reads, ect. stats gives a lengthy text report to be used with Multiqc. idxstats gives a table of read counts per sequence by chromosome, length, mapped_reads and unmapped_reads to give a distribution.\
\
To run the script on all of the BAM files, I used the following loop:

```{bash}
for bam in results/3_star/*Aligned.sortedByCoord.out.bam; do
sbatch scripts/samtools.sh $bam results/3_star results/3.5_samtools
done
```

# Salmon

Salmon is an RNA-seq quantification tool used to estimate transcript-level expression from RNA-seq reads. It assigns reads to transcripts and outputs expression measures such as TPM, estimated counts, and effective transcript lengths, which can later be summarized to the gene level for downstream analyses like differential expression, clustering, and exploratory data analysis. However, for Salmon to function correctly, it requires a reference transcriptome (FASTA of transcript sequences).

First, create a transcriptome index using gffread.

```{bash}
mkdir -p ref/transcripts
apptainer exec oras://community.wave.seqera.io/library/gffread:0.12.7--b08e770b84a4a126 \
  gffread ref/gtf/HelTub_1.0.gtf \
    -g ref/figshare/HelTub_1.0.fna \
    -w ref/transcripts/HelTub_1.0.transcripts.fa

```

Run the Salmon script using:

```{bash}
for bam in results/3_star/*Aligned.toTranscriptome.out.bam; do
sbatch scripts/salmon.sh $bam ref/transcripts/HelTub_1.0.transcripts.fa ref/gtf/HelTub_1.0.gtf results/4_salmon
done
```

# Final MultiQC

A comprehensive MultiQC can be ran using the raw fastq , trimmed fastq, star alignment, samtools results and salmon files to analyze the pipeline from start to finish and to determine if mapping rates and other statistics are optimal for downstream analysis.

To run the script, input:

```{bash}
sbatch scripts/complete_multiqc.sh results/1_raw_fastqc results/2_trimming results/3_star results/3.5_samtools results/4_salmon results/5_multiqc
```

Results:

| Sample Name | Error rate | Non-primary | Reads mapped | \% Proper pairs | Reads mapped | Aligned | Uniq aligned |
|:--------|:--------|:--------|:--------|:--------|:--------|:--------|:--------|
| SRR1555697  | 0.00%      | 31.0M       | 46.4M        | 100.0%          | 77.4M        | 91.7%   | 52.6%        |
| SRR1555695  | 0.00%      | 30.1M       | 43.7M        | 100.0%          | 73.7M        | 91.0%   | 51.1%        |
| SRR1555696  | 0.00%      | 30.8M       | 42.1M        | 100.0%          | 72.9M        | 92.2%   | 52.3%        |
| SRR1555736  | 0.00%      | 29.6M       | 41.0M        | 100.0%          | 70.6M        | 81.9%   | 46.5%        |
| SRR1555737  | 0.00%      | 28.3M       | 36.5M        | 100.0%          | 64.8M        | 73.2%   | 40.9%        |

## Interpretation

Mapping rates were really good! The lowest percentage that a read mapped was 73% which is good considering the plants are different cultivar, while the highest was 92% which is exceptionally good. Based on these results, the parameters seem correct for the scripts and analysis in R can take place.

# Prepping for R Analysis:

To transfer the necessary files we need from VScode to R studio, we need to first make a directory to that we can pull the files out of:

```         
mkdir results/6_R_analysis
```

In this file we will need:

1.  The quant files from salmon

2.  "samples.txt" which is the same as "runs.txt"

3.  A tx2gene.tsv file, which we will create last before moving to R studio.

    1.  Copy the quant files from results/4_salmon to results/6_R_analysis

        ```         
        cp results/4_salmon/SRR1555695/quant.sf results/6_R_analysis/SRR1555695_quant.sf 
        cp results/4_salmon/SRR1555696/quant.sf results/6_R_analysis/SRR1555696_quant.sf 
        cp results/4_salmon/SRR1555697/quant.sf results/6_R_analysis/SRR1555697_quant.sf 
        cp results/4_salmon/SRR1555736/quant.sf results/6_R_analysis/SRR1555636_quant.sf 
        cp results/4_salmon/SRR1555737/quant.sf results/6_R_analysis/SRR1555637_quant.sf 
        ```

    2.  Move and rename `runs.txt` to `samples.txt`

        ```         
        mv docs/runs.txt results/6_R_analysis/samples.txt
        ```

    3.  Create a file with gene-to-transfer attributes (below)

## Creating a file with gene-to-transfer attributes

Create a tsv file from the gtf file with column names

```{bash}
cut -f 9 ref/gtf/HelTub_1.0.gtf > gene2transferdraft.txt
```

`{"Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"} transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"; transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"; transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"; transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"; transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"; transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"; transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"; transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144"; transcript_id "Htub.h1tg000837l.g144.t1"; gene_id "Htub.h1tg000837l.g144";`

`transcript_id` `"Htub.h1tg000837l.g144.t1";` `gene_id` and `"Htub.h1tg000837l.g144";`are all seperated by spaces. We can change the spaces into tabs using `sed` and convert that into a tsv file.

### sed

The `sed` command can be used to edit files, and acts as the "find and replace tool" of bash.

sed uses the format:

`sed s/old text/new text/how often` or how I see it:

`sed s/find/replace with/replace all?`

The first function I will use is `s/ /\t/g`. This says "find a single space character ' ', replace with a tab character (\\t is the escape character for tab), and do so globally or our"replace all" option.

```{bash}
sed 's/ /\t/g' gene2transferdraft.txt > gene2transferdraft.tsv
```

Run a pipe to separate out the transcript_id and the gene_id, remove the quotation marks and semicolons and keep only unique values.

```{bash}
cut -f 2,4 gene2transferdraft.tsv | sed 's/"//g'| sed 's/;//g'| uniq > gene2transferdraft2.tsv
```

sed here is replacing he quotation marks (") with nothing (//) and doing so globally. The same is done for semicolon (;).

```{bash}
(echo -e "TXNAME\tGENEID"
  cat gene2transferdraft2.tsv) > tx2gene.tsv
```

echo -e:\
the -e option with echo prints escape characters, starting with `\` the `\t` option indicates to print tab, so `TXNAME\tGENEID` becomes `TXNAME     GENEID`seperated by a tab, which seperates into different columns in the tsv file.

check the file output

```{bash}
head /fs/ess/PAS2880/users/bateman139/project/ref/gtf/tx2gene.tsv
```

Move to ref dir and remove draft files

```{bash}
mv tx2gene.tsv results/6_R_analysis
rm *draft*
```

# 

# Continue To R Studio For Analysis
