#!/bin/bash
#SBATCH --partition=pibu_el8
#SBATCH --job-name=fastp_trim
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/FastP/logs/fastp_trim_%A_%a.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/FastP/logs/fastp_trim_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --mem=4000
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch


# DIRECTORIES (HOST)


WORKDIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood
RAW=$WORKDIR/raw_data
OUTDIR=$WORKDIR/FastP
LOGDIR=$OUTDIR/logs
CONTAINER=/containers/apptainer/fastp_0.24.1.sif

mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"


# DETECT R1 FILES

FILES=(${RAW}/*_1.fastq.gz)
N=${#FILES[@]}

if [ "$N" -eq 0 ]; then
    echo "ERROR: No R1 files found in $RAW"
    exit 1
fi


# RELAUNCH AS ARRAY

if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Submitting array with $N tasks..."
    sbatch --array=0-$(($N - 1)) "$0"
    exit 0
fi


# ASSIGN FILES FOR THIS TASK

R1=${FILES[$SLURM_ARRAY_TASK_ID]}
R2=${R1/_1.fastq.gz/_2.fastq.gz}
SAMPLE=$(basename "$R1" _1.fastq.gz)

OUT1=$OUTDIR/${SAMPLE}_1.clean.fastq.gz
OUT2=$OUTDIR/${SAMPLE}_2.clean.fastq.gz
JSON=$LOGDIR/${SAMPLE}.json
HTML=$LOGDIR/${SAMPLE}.html

echo "=== Processing sample: $SAMPLE ==="


# RUN FASTP inside APPTAINER

apptainer exec \
    --bind /data:/data \
    "$CONTAINER" fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$OUT1" \
        -O "$OUT2" \
        -j "$JSON" \
        -h "$HTML" \
        --detect_adapter_for_pe \
        -w $SLURM_CPUS_PER_TASK

echo "=== Finished sample: $SAMPLE ==="

