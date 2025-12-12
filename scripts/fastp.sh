#!/bin/bash
#SBATCH --account=PAS2880
#SBATCH --job-name=fastp_trim
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=8G

set -euo pipefail

# Set date so all outputs from this run are consistent (MMDD format)


run_date=$(date +%m%d)

# Positional arguments:
R1="$1"
R2="$2"
TRIM_DIR="$3"
FASTP_REPORT_DIR=$TRIM_DIR/fastp_report

#Logging
echo "Starting trimming script on $R1 and $R2"
echo "Saving trimmed files to $TRIM_DIR"
echo "----------------------------------------"

# Create directories automatically
mkdir -p "$TRIM_DIR" "$FASTP_REPORT_DIR"

    # Extract sample name for file naming
    # Example: data/fastq/SRR1555696_1.fastq.gz -> SRR1555696
    sample=${R1##*/}
    sample=${sample%_1.fastq.gz}

#Run Fastp
    apptainer exec \
      oras://community.wave.seqera.io/library/fastp:1.0.1--a5a7772c43b5ebcb \
      fastp \
      -i "$R1" \
      -I "$R2" \
      -o "${TRIM_DIR}/${sample}_${run_date}_1.trimmed.fastq.gz" \
      -O "${TRIM_DIR}/${sample}_${run_date}_2.trimmed.fastq.gz" \
      --cut_tail \
      --cut_tail_window_size 4 \
      --cut_tail_mean_quality 20 \
      -n 5 \
      --length_required 30 \
      --detect_adapter_for_pe \
      -h "${FASTP_REPORT_DIR}/${sample}_${run_date}.html" \
      -j "${FASTP_REPORT_DIR}/${sample}_${run_date}.json" \
      --thread 8
    # -i/-I: input reads
    # -o/-O: trimmed outputs (date-stamped)
    # --cut_tail* : 3' quality trimming with 4-bp window at Q20
    # -n 5: discard reads with >5 Ns
    # -l 50: discard reads shorter than 50 bp after trimming
    # --detect_adapter_for_pe: safe adapter detection
    # -h/-j: HTML/JSON reports (date-stamped)
    # --thread 8: use 8 threads


#Final Logging 
echo "Trimming on $R1 and $R2 complete"
