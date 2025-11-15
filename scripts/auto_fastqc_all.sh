#!/bin/bash
#SBATCH --partition=pshort_el8
#SBATCH --account=gmaringar+
#SBATCH --job-name=fastqc_array
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/fastqc_%A_%a.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/fastqc_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=1000
#SBATCH --cpus-per-task=1

# ============================
# AUTO-DETECT NUMBER OF FILES
# ============================

WORKDIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood
RAW=${WORKDIR}/raw_data
OUT=${WORKDIR}/qc/fastqc
mkdir -p ${OUT}

# Count number of FASTQ files
NFILES=$(ls ${RAW}/*.fastq.gz | wc -l)
echo "Detected $NFILES FASTQ files."

# If this script is not yet running as an array job, relaunch itself as one
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
  echo "Submitting array job with $NFILES tasks..."
  sbatch --array=0-$(($NFILES - 1)) $0
  exit 0
fi

# ============================
# MAIN EXECUTION
# ============================

module purge

CONTAINER=/containers/apptainer/fastqc-0.12.1.sif
FILES=(${RAW}/*.fastq.gz)
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

echo "Running FastQC on $FILE"
apptainer exec --bind /data/ ${CONTAINER} fastqc "$FILE" -o ${OUT}

echo "FastQC completed for $FILE"

