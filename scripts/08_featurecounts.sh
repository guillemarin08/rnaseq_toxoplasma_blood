#!/bin/bash

#SBATCH --job-name=featureCounts
#SBATCH --cpus-per-task=4
#SBATCH --mem=4000M
#SBATCH --time=2:00:00
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/featurecounts_%j.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/featurecounts_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch


# Paths

WORKDIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood
CONTAINER=/containers/apptainer/subread_2.0.6.sif

ANNOTATION=${WORKDIR}/genome/Mus_musculus.GRCm39.115.gtf
BAM_DIR=${WORKDIR}/mapped
OUT_DIR=${WORKDIR}/counts

mkdir -p $OUT_DIR


# Run featureCounts (single job for ALL BAM files)

apptainer exec --bind /data $CONTAINER \
    featureCounts \
        -T $SLURM_CPUS_PER_TASK \        # Number of threads (matched to SLURM settings)
        -a $ANNOTATION \                 # GTF annotation file (gene coordinates)
        -o ${OUT_DIR}/counts.txt \       # Output file containing gene counts
        -s 2 \                           # Strandedness: 2 = reversely stranded (Illumina TruSeq / RF)
        -p \                             # Paired-end mode: count fragments instead of individual reads
        ${BAM_DIR}/*_sorted.bam          # All sorted BAM files from HISAT2 mapping

echo "Finished FeatureCounts."

