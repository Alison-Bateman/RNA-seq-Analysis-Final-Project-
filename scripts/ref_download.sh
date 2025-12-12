#!/bin/bash
#SBATCH --account=PAS2880

set -euo pipefail

#Set Positional Parameters: 
LINK=$1
OUTDIR=$2
FILE_NAME=$3

#Initial logging: 
mkdir -p $OUTDIR
echo 
echo "Downloading $FILE_NAME from $LINK"
echo "Outdir: $OUTDIR"
wget \
  --content-disposition \
  --trust-server-names \
  --user-agent="Mozilla/5.0" \
  $LINK \
  -O $OUTDIR/$FILE_NAME
ls -lh $OUTDIR/$FILE_NAME

#final logging 
echo "---------------------------------------"
echo "download complete"