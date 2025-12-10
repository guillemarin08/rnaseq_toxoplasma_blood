#!/bin/bash
#SBATCH --job-name=hisat2_index
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/genome/hisat2_index/index_%j.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/genome/hisat2_index/index_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch


# Paths

CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif

GENOME_DIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/genome
FASTA_GZ=${GENOME_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
FASTA=${GENOME_DIR}/Mus_musculus.GRCm39.dna.primary_assembly.fa

INDEX_DIR=${GENOME_DIR}/hisat2_index

mkdir -p ${INDEX_DIR}


# Uncompress FASTA if needed

if [ ! -f "$FASTA" ]; then
    echo "=== Uncompressing FASTA ==="
    gunzip -c "$FASTA_GZ" > "$FASTA"
else
    echo "=== FASTA already uncompressed ==="
fi


# Build HISAT2 index

echo "=== Building HISAT2 index ==="

apptainer exec --bind /data "$CONTAINER" \
    hisat2-build -p $SLURM_CPUS_PER_TASK \
    "$FASTA" "${INDEX_DIR}/hisat2_index"

echo "=== Indexing complete ==="

