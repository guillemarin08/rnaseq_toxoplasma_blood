#!/bin/bash
#SBATCH --job-name=merge_counts
#SBATCH --partition=pibu_el8
#SBATCH --output=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/merge_counts_%j.out
#SBATCH --error=/data/users/gmaringarcia/rnaseq_toxoplasma_blood/logs/merge_counts_%j.err
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

COUNTS_DIR="/data/users/gmaringarcia/rnaseq_toxoplasma_blood/counts"
OUTPUT="$COUNTS_DIR/all_samples_counts.txt"

cd $COUNTS_DIR

FILES=( $(ls *.counts.txt) )

echo "Found ${#FILES[@]} files."

# ----- Extract Gene IDs from the first file (skip comments, skip second header) -----

FIRST=${FILES[0]}

grep -v "^#" "$FIRST" | sed '1d' | cut -f1 > gene_ids.tmp

# Start matrix
cp gene_ids.tmp matrix.tmp

# ----- Loop through each sample -----
for FILE in "${FILES[@]}"; do

    SAMPLE=$(basename "$FILE" .counts.txt)

    echo "Processing $SAMPLE"

    # Skip comment lines (#...) AND skip the header line with BAM paths
    grep -v "^#" "$FILE" | sed '1d' | awk '{print $NF}' > ${SAMPLE}.tmp

    paste matrix.tmp ${SAMPLE}.tmp > matrix2.tmp
    mv matrix2.tmp matrix.tmp

done

# ----- Create clean header -----
HEADER="Geneid"
for FILE in "${FILES[@]}"; do
    SAMPLE=$(basename "$FILE" .counts.txt)
    HEADER="${HEADER}\t${SAMPLE}"
done

echo -e "$HEADER" > "$OUTPUT"
cat matrix.tmp >> "$OUTPUT"

rm *.tmp

echo "Matrix created at: $OUTPUT"
