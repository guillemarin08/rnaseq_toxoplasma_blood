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
outdir <- file.path(project_dir, "deseq2/A_WT_case_vs_control")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#==============================
# LOAD DATA
#==============================

counts <- read.table(counts_file, header = TRUE, row.names = 1)
meta <- read.table(metadata_file, header = TRUE, sep = "\t")

#==============================
# FILTER WT ONLY
#==============================

meta_wt <- meta[meta$genotype == "WT", ]
counts_wt <- counts[, meta_wt$sample]

# Make sure infection has correct reference level
meta_wt$infection <- factor(meta_wt$infection, levels = c("Control", "Case"))

#==============================
# CREATE DESEQ2 OBJECT
#==============================

dds_wt <- DESeqDataSetFromMatrix(
  countData = counts_wt,
  colData = meta_wt,
  design = ~ infection
)

dds_wt <- DESeq(dds_wt)

#==============================
# RESULTS: Case vs Control
#==============================

res_A <- results(dds_wt, contrast = c("infection", "Case", "Control"))
res_A <- res_A[order(res_A$padj), ]

write.csv(as.data.frame(res_A),
          file.path(outdir, "DEG_A_WT_case_vs_control.csv"))

#==============================
# VST + QC PLOTS
#==============================

vsd <- vst(dds_wt)

# PCA plot
pca_data <- plotPCA(vsd, intgroup = "infection", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = infection)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_bw()

ggsave(file.path(outdir, "PCA_A_WT_case_vs_control.png"), p,
       width = 6, height = 5)

# MA plot
png(file.path(outdir, "MAplot_A_WT_case_vs_control.png"), width = 1000, height = 800)
plotMA(res_A, ylim = c(-5, 5))
dev.off()

# Heatmap top 50
top50 <- head(rownames(res_A), 50)
mat <- assay(vsd)[top50, ]
mat <- mat - rowMeans(mat)

anno <- data.frame(condition = meta_wt$infection)
rownames(anno) <- meta_wt$sample

pheatmap(mat,
         annotation_col = anno,
         filename = file.path(outdir, "heatmap_A_top50.png"))

# Save DESeq2 object
saveRDS(dds_wt, file.path(outdir, "dds_wt_A.rds"))
