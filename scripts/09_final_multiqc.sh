#!/bin/bash
#SBATCH --job-name=multiqc_final
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/multiqc/final_multiqc.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/multiqc/final_multiqc.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch


# Paths

WORKDIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood
CONTAINER=/containers/apptainer/multiqc-1.19.sif

FASTQC_TRIM=$WORKDIR/FastP/FastP_QC
FASTP_LOGS=$WORKDIR/FastP/logs
MAPPING=$WORKDIR/mapped
FEATURECOUNTS=$WORKDIR/counts

OUTDIR=$WORKDIR/multiqc/final
mkdir -p $OUTDIR


# Run MultiQC

apptainer exec --bind /data $CONTAINER multiqc \
    $FASTQC_TRIM \
    $FASTP_LOGS \
    $MAPPING \
    $FEATURECOUNTS \
    -o $OUTDIR \
    --force

echo "Final MultiQC report created in: $OUTDIR"

