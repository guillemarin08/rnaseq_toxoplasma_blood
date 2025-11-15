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
outdir <- file.path(project_dir, "deseq2/C_WT_vs_DKO_control")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

#==============================
# LOAD DATA
#==============================

counts <- read.table(counts_file, header = TRUE, row.names = 1)
meta <- read.table(metadata_file, header = TRUE, sep = "\t")

#==============================
# FILTER: ONLY CONTROL SAMPLES FROM BOTH GENOTYPES
#==============================

meta_sub <- meta[meta$infection == "Control", ]
counts_sub <- counts[, meta_sub$sample]

# Set factor levels (WT must be the reference!)
meta_sub$genotype <- factor(meta_sub$genotype, levels = c("WT", "DKO"))

#==============================
# CREATE DESEQ2 OBJECT
#==============================

dds_c <- DESeqDataSetFromMatrix(
  countData = counts_sub,
  colData = meta_sub,
  design = ~ genotype
)

dds_c <- DESeq(dds_c)

#==============================
# RESULTS: DKO vs WT (baseline)
#==============================

res_C <- results(dds_c, contrast = c("genotype", "DKO", "WT"))
res_C <- res_C[order(res_C$padj), ]

write.csv(as.data.frame(res_C),
          file.path(outdir, "DEG_C_WT_vs_DKO_control.csv"))

#==============================
# VST + QC PLOTS
#==============================

vsd <- vst(dds_c)

# PCA
pca_data <- plotPCA(vsd, intgroup = "genotype", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

p <- ggplot(pca_data, aes(PC1, PC2, color = genotype)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  theme_bw()

ggsave(file.path(outdir, "PCA_C_WT_vs_DKO_control.png"), p,
       width = 6, height = 5)

# MA plot
png(file.path(outdir, "MAplot_C_WT_vs_DKO_control.png"), width = 1000, height = 800)
plotMA(res_C, ylim = c(-5, 5))
dev.off()

# Heatmap top 50
top50 <- head(rownames(res_C), 50)
mat <- assay(vsd)[top50, ]
mat <- mat - rowMeans(mat)

anno <- data.frame(genotype = meta_sub$genotype)
rownames(anno) <- meta_sub$sample

pheatmap(mat,
         annotation_col = anno,
         filename = file.path(outdir, "heatmap_C_top50.png"))

# Save DESeq2 object
saveRDS(dds_c, file.path(outdir, "dds_c_C.rds"))
