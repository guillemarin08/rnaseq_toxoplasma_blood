#!/usr/bin/env Rscript

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggplot2)

#==============================
# Paths
#==============================

project_dir <- "/data/users/gmaringarcia/rnaseq_toxoplasma_blood"
counts_file <- file.path(project_dir, "counts/all_samples_counts.txt")
metadata_file <- file.path(project_dir, "metadata/sample_metadata.tsv")
outdir <- file.path(project_dir, "deseq2/B_DKO_case_vs_control")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#==============================
# LOAD DATA
#==============================

counts <- read.table(counts_file, header = TRUE, row.names = 1)
meta <- read.table(metadata_file, header = TRUE, sep = "\t")

#==============================
# FILTER DKO ONLY
#==============================

meta_dko <- meta[meta$genotype == "DKO", ]
counts_dko <- counts[, meta_dko$sample]

# Set factor levels explicitly
meta_dko$infection <- factor(meta_dko$infection, levels = c("Control", "Case"))

#==============================
# CREATE DESEQ2 OBJECT
#==============================

dds_dko <- DESeqDataSetFromMatrix(
  countData = counts_dko,
  colData = meta_dko,
  design = ~ infection
)

dds_dko <- DESeq(dds_dko)

#==============================
# RESULTS: Case vs Control
#==============================

res_B <- results(dds_dko, contrast = c("infection", "Case", "Control"))
res_B <- res_B[order(res_B$padj), ]

write.csv(as.data.frame(res_B),
          file.path(outdir, "DEG_B_DKO_case_vs_control.csv"))

#==============================
# VST + QC PLOTS
#==============================

vsd <- vst(dds_dko)

# PCA
pca_data <- plotPCA(vsd, intgroup = "infection", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = infection)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_bw()

ggsave(file.path(outdir, "PCA_B_DKO_case_vs_control.png"), p,
       width = 6, height = 5)

# MA plot
png(file.path(outdir, "MAplot_B_DKO_case_vs_control.png"), width = 1000, height = 800)
plotMA(res_B, ylim = c(-5, 5))
dev.off()

# Heatmap top 50
top50 <- head(rownames(res_B), 50)
mat <- assay(vsd)[top50, ]
mat <- mat - rowMeans(mat)

anno <- data.frame(infection = meta_dko$infection)
rownames(anno) <- meta_dko$sample

pheatmap(mat,
         annotation_col = anno,
         filename = file.path(outdir, "heatmap_B_top50.png"))

# Save DESeq2 object
saveRDS(dds_dko, file.path(outdir, "dds_dko_B.rds"))
