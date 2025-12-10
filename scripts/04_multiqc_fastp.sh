#!/usr/bin/env bash

#SBATCH --job-name=multiqc_fastp_clean
#SBATCH --cpus-per-task=1
#SBATCH --mem=1000
#SBATCH --time=0:30:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/FastP/FastP_QC/MultiQC_clean/output_%x_%j.o
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/FastP/FastP_QC/MultiQC_clean/error_%x_%j.e


# Paths

WORKDIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood
CONTAINER=/containers/apptainer/multiqc-1.19.sif

INPUT_DIR=$WORKDIR/FastP/FastP_QC
OUTPUT_DIR=$WORKDIR/FastP/FastP_QC/MultiQC_clean


# Create output directories

mkdir -p "$OUTPUT_DIR"

echo "Running MultiQC on: $INPUT_DIR"


# Run MultiQC inside Apptainer

apptainer exec --bind $WORKDIR \
    "$CONTAINER" multiqc "$INPUT_DIR" -o "$OUTPUT_DIR" --force

echo "MultiQC report completed!"

