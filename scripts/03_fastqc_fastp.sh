#!/bin/bash
#SBATCH --partition=pshort_el8
#SBATCH --job-name=fastqc_fastp
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/FastP/FastP_QC/fastqc_%A_%a.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/FastP/FastP_QC/fastqc_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=1000
#SBATCH --cpus-per-task=1


# Directories

WORKDIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood
FASTP_DIR=$WORKDIR/FastP
OUTDIR=$FASTP_DIR/FastP_QC

mkdir -p "$OUTDIR"


# Detect cleaned FASTQ files (*.clean.fastq.gz)

FILES=(${FASTP_DIR}/*clean.fastq.gz)
NFILES=${#FILES[@]}

if [ "$NFILES" -eq 0 ]; then
  echo "ERROR: No cleaned FASTQ files found in $FASTP_DIR"
  exit 1
fi

echo "Detected $NFILES cleaned FASTQ files."


# Relaunch as array job if needed

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  echo "Submitting array job with $NFILES tasks..."
  sbatch --array=0-$(($NFILES - 1)) "$0"
  exit 0
fi


# Main execution

module purge
CONTAINER=/containers/apptainer/fastqc-0.12.1.sif

FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

echo "Running FastQC on: $FILE"

apptainer exec --bind /data/ $CONTAINER fastqc "$FILE" -o "$OUTDIR"

echo "FastQC completed for: $FILE"

