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

# Paths
COUNTS_DIR="/data/users/gmaringarcia/rnaseq_toxoplasma_blood/counts"
OUTPUT="${COUNTS_DIR}/all_samples_counts.txt"

cd $COUNTS_DIR

# Detect .counts.txt files
FILES=( $(ls *.counts.txt) )

echo "Detected ${#FILES[@]} files:"
printf '%s\n' "${FILES[@]}"

# Use the first file to extract gene IDs
FIRST=${FILES[0]}

echo "Using $FIRST as reference."

# Extract Geneid column
cut -f1 "$FIRST" > gene_column.tmp

# For each file, extract counts and add to matrix
for FILE in "${FILES[@]}"; do
    SAMPLE=$(basename "$FILE" .counts.txt)
    echo "Processing $SAMPLE"

    # Extract last column (= counts)
    cut -f7 "$FILE" > ${SAMPLE}.tmp

    # Paste into matrix
    paste gene_column.tmp ${SAMPLE}.tmp > temp && mv temp gene_column.tmp

    # Rename column inside header later
done

# Build header
HEADER="Geneid"
for FILE in "${FILES[@]}"; do
    SAMPLE=$(basename "$FILE" .counts.txt)
    HEADER="${HEADER}\t${SAMPLE}"
done

# Write output
echo -e "$HEADER" > "$OUTPUT"
cat gene_column.tmp >> "$OUTPUT"

# Clean temporary files
rm *.tmp

echo "Matrix created at:"
echo "$OUTPUT"
