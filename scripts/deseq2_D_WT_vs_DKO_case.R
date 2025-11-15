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
outdir <- file.path(project_dir, "deseq2/D_WT_vs_DKO_case")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#==============================
# LOAD DATA
#==============================

counts <- read.table(counts_file, header = TRUE, row.names = 1)
meta <- read.table(metadata_file, header = TRUE, sep = "\t")

#==============================
# FILTER: ONLY CASE SAMPLES FROM BOTH GENOTYPES
#==============================

meta_sub <- meta[meta$infection == "Case", ]
counts_sub <- counts[, meta_sub$sample]

# Make sure WT is reference level
meta_sub$genotype <- factor(meta_sub$genotype, levels = c("WT", "DKO"))

#==============================
# CREATE DESEQ2 OBJECT
#==============================

dds_d <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData = meta_sub,
  design = ~ genotype
)

dds_d <- DESeq(dds_d)

#==============================
# RESULTS: DKO vs WT (Case condition)
#==============================

res_D <- results(dds_d, contrast = c("genotype", "DKO", "WT"))
res_D <- res_D[order(res_D$padj), ]

write.csv(as.data.frame(res_D),
          file.path(outdir, "DEG_D_WT_vs_DKO_case.csv"))

#==============================
# VST + QC PLOTS
#==============================

vsd <- vst(dds_d)

# PCA
pca_data <- plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = genotype)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_bw()

ggsave(file.path(outdir, "PCA_D_WT_vs_DKO_case.png"), p,
       width = 6, height = 5)

# MA plot
png(file.path(outdir, "MAplot_D_WT_vs_DKO_case.png"), width = 1000, height = 800)
plotMA(res_D, ylim = c(-5, 5))
dev.off()

# Heatmap top 50
top50 <- head(rownames(res_D), 50)
mat <- assay(vsd)[top50, ]
mat <- mat - rowMeans(mat)

anno <- data.frame(genotype = meta_sub$genotype)
rownames(anno) <- meta_sub$sample

pheatmap(mat,
         annotation_col = anno,
         filename = file.path(outdir, "heatmap_D_top50.png"))

# Save DESeq2 object
saveRDS(dds_d, file.path(outdir, "dds_d_D.rds"))
