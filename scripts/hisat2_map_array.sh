#!/bin/bash
#SBATCH --job-name=hisat2_map
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/hisat2_map_%A_%a.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/hisat2_map_%A_%a.err
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=0-14
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch

# -------------------------------------------------------------
#   No "module load apptainer" — apptainer is globally available
# -------------------------------------------------------------

CONTAINER="/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif"

RAW_DIR="/data/users/gmaringarcia/rnaseq_toxoplasma_blood/trimmed"
INDEX_DIR="/data/users/gmaringarcia/rnaseq_toxoplasma_blood/genome/hisat2_index"
OUT_DIR="/data/users/gmaringarcia/rnaseq_toxoplasma_blood/mapped"

mkdir -p $OUT_DIR

SAMPLES=($(ls ${RAW_DIR}/*_1.trimmed.fastq.gz))
R1=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
R2=${R1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}

BASENAME=$(basename $R1 _1.trimmed.fastq.gz)

SAM=${OUT_DIR}/${BASENAME}.sam
BAM=${OUT_DIR}/${BASENAME}.bam
SORTED=${OUT_DIR}/${BASENAME}.sorted.bam

# 1. Align
apptainer exec --bind /data ${CONTAINER} \
    hisat2 -x ${INDEX_DIR}/genome_index \
    -1 $R1 -2 $R2 \
    --rna-strandness RF \
    -p $SLURM_CPUS_PER_TASK \
    -S $SAM

# 2. SAM → BAM
apptainer exec --bind /data ${CONTAINER} \
    samtools view -hbS $SAM > $BAM
rm $SAM

# 3. Sort
apptainer exec --bind /data ${CONTAINER} \
    samtools sort -@ $SLURM_CPUS_PER_TASK \
    -o $SORTED \
    -T ${OUT_DIR}/${BASENAME}.tmp \
    $BAM

rm $BAM

# 4. Index
apptainer exec --bind /data ${CONTAINER} \
    samtools index $SORTED

echo "Done: $BASENAME"

