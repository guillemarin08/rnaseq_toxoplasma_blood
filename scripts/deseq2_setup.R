library(DESeq2)
library(tidyverse)

# Load count matrix
counts <- read.table(
  "/data/users/gmaringarcia/rnaseq_toxoplasma_blood/counts/all_counts.tsv",
  header = TRUE, row.names = 1
)

# Load metadata
meta <- read.table(
  "/data/users/gmaringarcia/rnaseq_toxoplasma_blood/metadata/sample_metadata.tsv",
  header = TRUE, sep = "\t"
)

# Ensure same sample order
all(colnames(counts) == meta$sample)
# Should print TRUE

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ genotype + infection
)

saveRDS(dds, "dds_raw.rds")

