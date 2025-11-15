#!/bin/bash
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastp_trim
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/fastp_trim_%A_%a.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/fastp_trim_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=4000
#SBATCH --cpus-per-task=4

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch

#############################################
# Load required modules
#############################################
module purge
module load fastp/0.23.4-GCC-10.3.0   # <-- ESTE ES EL MÓDULO DE TU CLÚSTER

#############################################
# Define directories
#############################################
WORKDIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood
RAW=${WORKDIR}/raw_data
TRIM=${WORKDIR}/trimmed
LOGS=${TRIM}/fastp_logs

mkdir -p "$TRIM" "$LOGS"

#############################################
# Detect R1 files
#############################################
FILES=(${RAW}/*_1.fastq.gz)
N=${#FILES[@]}

if [ "$N" -eq 0 ]; then
  echo "ERROR: No R1 FASTQ files found in $RAW using pattern *_1.fastq.gz"
  exit 1
fi

#############################################
# Relaunch script as SLURM array if needed
#############################################
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Submitting SLURM array job with $N tasks..."
    sbatch --array=0-$(($N - 1)) "$0"
    exit 0
fi

#############################################
# Assign R1 and R2 for this task
#############################################
R1=${FILES[$SLURM_ARRAY_TASK_ID]}
R2=${R1/_1.fastq.gz/_2.fastq.gz}

if [ ! -f "$R2" ]; then
  echo "ERROR: Matching R2 file not found for $R1"
  exit 1
fi

BASENAME=$(basename "$R1" _1.fastq.gz)

OUT1=${TRIM}/${BASENAME}_1.trimmed.fastq.gz
OUT2=${TRIM}/${BASENAME}_2.trimmed.fastq.gz
JSON=${LOGS}/${BASENAME}.json
HTML=${LOGS}/${BASENAME}.html

echo "=== Running fastp for sample: $BASENAME ==="

#############################################
# Run fastp
#############################################
fastp \
    -i "$R1" -I "$R2" \
    -o "$OUT1" -O "$OUT2" \
    --detect_adapter_for_pe \
    --cut_front \
    --cut_tail \
    --cut_mean_quality 20 \
    --length_required 35 \
    --thread "$SLURM_CPUS_PER_TASK" \
    --html "$HTML" \
    --json "$JSON"

echo "=== Finished sample: $BASENAME ==="

