#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/hisat2_index_%j.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/hisat2_index_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch
#SBATCH --mail-type=BEGIN,END,FAIL

# Load container (no module needed)
CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

# Paths
GENOME_DIR="/data/users/gmaringarcia/rnaseq_toxoplasma_blood/genome"
FASTA="${GENOME_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
INDEX_DIR="${GENOME_DIR}/hisat2_index"

mkdir -p ${INDEX_DIR}

echo "=== Uncompressing FASTA ==="
gunzip -c $FASTA > ${GENOME_DIR}/genome.fa

echo "=== Building HISAT2 index ==="
apptainer exec --bind /data/ $CONTAINER \
    hisat2-build -p 8 ${GENOME_DIR}/genome.fa ${INDEX_DIR}/genome_index

echo "=== Done ==="

