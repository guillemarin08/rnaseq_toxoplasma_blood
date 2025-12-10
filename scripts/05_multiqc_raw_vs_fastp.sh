#!/bin/bash
#SBATCH --job-name=multiqc_raw_vs_fastp
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/multiqc/multiqc_raw_vs_fastp/multiqc_%j.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/multiqc/multiqc_raw_vs_fastp/multiqc_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch


# Paths

WORKDIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood

RAW_QC=$WORKDIR/qc/fastqc
FASTP_QC=$WORKDIR/FastP/FastP_QC
FASTP_LOGS=$WORKDIR/FastP/logs

OUTDIR=$WORKDIR/multiqc/multiqc_raw_vs_fastp
mkdir -p $OUTDIR

# MultiQC container
CONTAINER=/containers/apptainer/multiqc-1.19.sif


# Run MultiQC using Apptainer

apptainer exec --bind /data ${CONTAINER} multiqc \
    $RAW_QC \
    $FASTP_QC \
    $FASTP_LOGS \
    -o $OUTDIR \
    --force

echo "MultiQC report created in: $OUTDIR"

