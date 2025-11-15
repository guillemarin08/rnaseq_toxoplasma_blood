#!/bin/bash
#SBATCH --job-name=multiqc_trim
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/multiqc_trim_%j.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/multiqc_trim_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch

# No need to load apptainer module – not present in pibu_el8
# Path to the MultiQC container
CONTAINER=/containers/apptainer/multiqc-1.19.sif

# Run MultiQC (fastp + fastqc)
apptainer exec --bind /data ${CONTAINER} multiqc \
    /data/users/gmaringarcia/rnaseq_toxoplasma_blood/trimmed \
    /data/users/gmaringarcia/rnaseq_toxoplasma_blood/trimmed/fastp_logs \
    -o /data/users/gmaringarcia/rnaseq_toxoplasma_blood/qc/multiqc_trim

