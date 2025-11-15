#!/bin/bash
#SBATCH --partition=pshort_el8
#SBATCH --job-name=multiqc
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/multiqc_%j.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/multiqc_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=2000
#SBATCH --cpus-per-task=1

module purge
module load MultiQC/1.11-foss-2021a

WORKDIR=/data/users/gmaringarcia/rnaseq_toxoplasma_blood
FASTQC_DIR=${WORKDIR}/qc/fastqc
OUT_DIR=${WORKDIR}/qc/multiqc

mkdir -p ${OUT_DIR}

NFILES=$(ls ${FASTQC_DIR}/*.zip 2>/dev/null | wc -l)
if [ "$NFILES" -eq 0 ]; then
  echo "❌ No FastQC .zip files found in ${FASTQC_DIR}. Run FastQC first."
  exit 1
fi

echo "Detected $NFILES FastQC result files. Running MultiQC..."

multiqc ${FASTQC_DIR} -o ${OUT_DIR}

echo "✅ MultiQC completed. Report saved in:"
echo "   ${OUT_DIR}/multiqc_report.html"

