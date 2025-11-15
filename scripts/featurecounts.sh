#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/featurecounts_%A_%a.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/featurecounts_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-14
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch

###########################################
# Load Subread module (contains featureCounts)
###########################################
module load Subread/2.0.3-GCC-10.3.0

###########################################
# Paths
###########################################
BAM_DIR="/data/users/gmaringarcia/rnaseq_toxoplasma_blood/mapped"
OUT_DIR="/data/users/gmaringarcia/rnaseq_toxoplasma_blood/counts"
GTF="/data/users/gmaringarcia/rnaseq_toxoplasma_blood/genome/Mus_musculus.GRCm39.115.gtf"

mkdir -p $OUT_DIR

###########################################
# Get list of BAM files
###########################################
BAM_FILES=($(ls ${BAM_DIR}/*.sorted.bam))
BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID]}"

# Extract sample name
SAMPLE=$(basename "$BAM" .sorted.bam)

echo "Processing sample: $SAMPLE"
echo "Input BAM: $BAM"

###########################################
# Run featureCounts
###########################################
featureCounts \
    -T 4 \
    -p \
    -s 2 \
    -a "$GTF" \
    -o "${OUT_DIR}/${SAMPLE}.counts.txt" \
    "$BAM"

echo "Done: ${SAMPLE}"

