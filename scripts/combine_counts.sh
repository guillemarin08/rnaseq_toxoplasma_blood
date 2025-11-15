#!/bin/bash
#SBATCH --job-name=merge_counts
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/merge_counts_%j.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/merge_counts_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=guillermo.maringarcia@students.unibe.ch

COUNTS_DIR="/data/users/gmaringarcia/rnaseq_toxoplasma_blood/counts"
OUTPUT="$COUNTS_DIR/all_samples_counts.txt"

cd $COUNTS_DIR

FILES=( $(ls *.counts.txt) )

echo "Merging ${#FILES[@]} count files."

FIRST=${FILES[0]}

# extract gene IDs (skip comment lines)
grep -v "^#" "$FIRST" | cut -f1 > gene_ids.tmp

# init final matrix with gene IDs
cp gene_ids.tmp matrix.tmp

# loop through samples
for FILE in "${FILES[@]}"; do
    SAMPLE=$(basename "$FILE" .counts.txt)

    echo "Processing $SAMPLE"

    # skip comment lines, extract last column (counts)
    grep -v "^#" "$FILE" | awk '{print $NF}' > ${SAMPLE}.tmp

    paste matrix.tmp ${SAMPLE}.tmp > matrix2.tmp
    mv matrix2.tmp matrix.tmp
done

# write header
HEADER="Geneid"
for FILE in "${FILES[@]}"; do
    SAMPLE=$(basename "$FILE" .counts.txt)
    HEADER="${HEADER}\t${SAMPLE}"
done

echo -e "$HEADER" > "$OUTPUT"
cat matrix.tmp >> "$OUTPUT"

rm *.tmp

echo "Final matrix created at: $OUTPUT"
