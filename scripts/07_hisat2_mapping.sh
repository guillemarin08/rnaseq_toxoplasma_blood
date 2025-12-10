#!/bin/bash
#SBATCH --job-name=hisat2_map
#SBATCH --partition=pibu_el8
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/mapping_%A_%a.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/mapping_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch


# Paths

CONTAINER=/containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif
DATADIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood
RAW_DIR=$DATADIR/FastP
INDEX_DIR=$DATADIR/genome/hisat2_index
OUT_DIR=$DATADIR/mapped

mkdir -p $OUT_DIR


# Autodetect samples

R1_FILES=(${RAW_DIR}/*_1.clean.fastq.gz)
N=${#R1_FILES[@]}

# Auto-launch array if needed
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    sbatch --array=0-$(($N-1)) $0
    exit 0
fi

R1=${R1_FILES[$SLURM_ARRAY_TASK_ID]}
R2=${R1/_1.clean.fastq.gz/_2.clean.fastq.gz}

SAMPLE=$(basename "$R1" _1.clean.fastq.gz)

echo "=== Processing sample: $SAMPLE ==="


# Output paths

SCRATCH_SAM=$SCRATCH/${SAMPLE}.sam
SORTED_BAM=${OUT_DIR}/${SAMPLE}_sorted.bam
SUMMARY=${OUT_DIR}/${SAMPLE}_summary.txt


# HISAT2 mapping

apptainer exec --bind /data $CONTAINER \
    hisat2 -p $SLURM_CPUS_PER_TASK \
           -x ${INDEX_DIR}/hisat2_index \
           -1 $R1 -2 $R2 \
           --rna-strandness RF \
           --summary-file $SUMMARY \
           -S $SCRATCH_SAM

echo "Mapping finished for $SAMPLE."


# SAM → sorted BAM

apptainer exec --bind /data $CONTAINER \
    samtools sort -@ $SLURM_CPUS_PER_TASK \
                  -o $SORTED_BAM \
                  $SCRATCH_SAM

rm $SCRATCH_SAM

echo "Sorted BAM saved: $SORTED_BAM"


# BAM index

apptainer exec --bind /data $CONTAINER \
    samtools index $SORTED_BAM

echo "Index created."
echo "=== DONE: $SAMPLE ==="

