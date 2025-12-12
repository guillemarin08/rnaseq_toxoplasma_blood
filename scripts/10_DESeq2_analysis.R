#!/usr/bin/env Rscript

# ==============================================================================
# RNA-seq Analysis Workflow: Interaction Effect (Genotype x Infection)
# ==============================================================================
# Description:
#   This script performs differential expression analysis using DESeq2.
#   It focuses on the interaction between Genotype (DKO vs WT) and Infection.
#   Workflow includes: Quality Control (PCA), DE Analysis, Visualization (Volcano, 
#   Counts), and Functional Enrichment (GO).
#
# Inputs:
#   - FeatureCounts matrix (tab-separated)
#   - Metadata file (tab-separated)
#
# Outputs:
#   - CSV tables for DE results and GO enrichment
#   - High-resolution PNG images for PCA, Volcano, Gene Counts, and GO plots.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Load Required Packages
# ------------------------------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db) # Mouse annotation database

# ------------------------------------------------------------------------------
# 1. Define Directories and Files
# ------------------------------------------------------------------------------
base_dir <- "C:/Users/Win10/Desktop/Bioinformatica máster/RNA seq/mis resultados"

# Locate input files automatically based on patterns
count_file <- list.files(base_dir, pattern = "^counts_deseq", full.names = TRUE)
meta_file  <- list.files(base_dir, pattern = "^sample_metadata", full.names = TRUE)

# ------------------------------------------------------------------------------
# 2. Load Data
# ------------------------------------------------------------------------------
counts <- read.delim(count_file, row.names = 1, check.names = FALSE)
coldata <- read.delim(meta_file,  row.names = 1, check.names = FALSE)

# Pre-processing: Remove non-count columns (Chr, Start, Length, etc.) if present
if("Length" %in% colnames(counts) | "Chr" %in% colnames(counts)) {
  counts <- counts[ , -c(1:5)]
}

# ------------------------------------------------------------------------------
# 3. Data Integrity Checks
# ------------------------------------------------------------------------------
# Ensure columns in counts match rows in metadata
if (!all(colnames(counts) %in% rownames(coldata))) {
  stop("ERROR: Mismatch between count matrix columns and metadata rows.")
}
coldata <- coldata[colnames(counts), ] # Align order

# 4. Define Factors and Reference Levels
coldata$genotype  <- factor(coldata$genotype)
coldata$infection <- factor(coldata$infection)

# Set Reference Levels: WT and Control
coldata$infection <- relevel(coldata$infection, ref = "Control")
coldata$genotype  <- relevel(coldata$genotype,  ref = "WT")

# ------------------------------------------------------------------------------
# 5. Run DESeq2 (Interaction Design)
# ------------------------------------------------------------------------------
# Design formula: ~ genotype * infection
# This tests for the interaction effect (difference in response between genotypes)
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ genotype * infection)

dds <- DESeq(dds)

# ------------------------------------------------------------------------------
# 6. Exploratory Data Analysis (PCA)
# ------------------------------------------------------------------------------
vsd <- vst(dds, blind = TRUE) # Variance stabilizing transformation for QC
pcaData <- plotPCA(vsd, intgroup = c("genotype", "infection"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p_pca <- ggplot(pcaData, aes(PC1, PC2, color = infection, shape = genotype)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA Plot: Sample Clustering") +
  theme_bw()

print(p_pca)

# Save Plot
ggsave(file.path(base_dir, "EDA_PCA_Plot.png"), plot = p_pca, width = 8, height = 6)


# ------------------------------------------------------------------------------
# 7. Differential Expression Analysis (Interaction Term)
# ------------------------------------------------------------------------------
# Extract results for the interaction term (genotypeDKO.infectionCase)
res_interaction <- results(dds, name = grep("genotype.*infection.*Case", resultsNames(dds), value=TRUE), alpha = 0.05)
summary(res_interaction)

# Filter significant genes
sig_genes <- res_interaction[which(res_interaction$padj < 0.05), ]
up_genes <- sig_genes[which(sig_genes$log2FoldChange > 0), ]
down_genes <- sig_genes[which(sig_genes$log2FoldChange < 0), ]

# Annotation: Map Ensembl IDs to Gene Symbols for readability
res_interaction$symbol <- mapIds(org.Mm.eg.db,
                                 keys = rownames(res_interaction),
                                 column = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")

sig_genes$symbol <- mapIds(org.Mm.eg.db,
                           keys = rownames(sig_genes),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")

# ------------------------------------------------------------------------------
# 8. Visualization: Volcano Plot 
# ------------------------------------------------------------------------------
# Style: Blue for Down-regulated, Red for Up-regulated, Grey for NS.

genes_highlight <- c("Tlr11", "Cxcl10", "Tnfsf9") # Specific genes to label

# Define custom color palette
keyvals <- rep('grey75', nrow(res_interaction))
names(keyvals) <- rep('Not significant', nrow(res_interaction))

# Down-regulated (Blue)
down_idx <- which(res_interaction$padj < 0.05 & res_interaction$log2FoldChange < -1.0)
keyvals[down_idx] <- 'blue'
names(keyvals)[down_idx] <- 'Down-regulated'

# Up-regulated (Red)
up_idx <- which(res_interaction$padj < 0.05 & res_interaction$log2FoldChange > 1.0)
keyvals[up_idx] <- 'red'
names(keyvals)[up_idx] <- 'Up-regulated'

p_volcano <- EnhancedVolcano(res_interaction,
                             lab = res_interaction$symbol,
                             x = 'log2FoldChange',
                             y = 'padj',
                             selectLab = genes_highlight,
                             xlab = bquote(~Log[2]~ 'fold change'),
                             ylab = bquote(~-Log[10]~ 'adjusted p-value'),
                             title = 'Genotype x Infection Interaction Landscape',
                             subtitle = 'Genes showing significantly different responses in DKO compared to WT',
                             pCutoff = 0.05,
                             FCcutoff = 1.0,
                             pointSize = 3.0,
                             labSize = 10.0,
                             colCustom = keyvals,
                             colAlpha = 0.8,
                             legendPosition = 'right',
                             legendLabSize = 16,
                             legendIconSize = 6.0,
                             drawConnectors = TRUE,
                             widthConnectors = 0.5,
                             boxedLabels = TRUE,
                             gridlines.major = FALSE, 
                             gridlines.minor = FALSE
) +
  theme_bw(base_size = 22) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    plot.subtitle = element_text(size = 24),
    axis.title = element_text(size = 25),
    axis.text  = element_text(size = 19),
    legend.title = element_blank(),               
    legend.text  = element_text(size = 23)
  ) +
  guides(color = guide_legend(title = NULL))
  

print(p_volcano)

#Save volcano plot
ggsave(file.path(base_dir, "DE_Volcano_Plot.png"), p_volcano, width = 13, height = 11, dpi = 300)

# ------------------------------------------------------------------------------
# 9. Visualization: Normalized Counts for Selected Genes
# ------------------------------------------------------------------------------
# Investigate expression levels for Tlr11, Cxcl10 and Tnfsf9 to visualize the interaction.

for (gene_sym in genes_highlight) {
  
  # Find Ensembl ID
  gene_id <- rownames(res_interaction)[which(res_interaction$symbol == gene_sym)]
  if(length(gene_id) > 0) gene_id <- gene_id[1] else next
  
  # Get normalized counts
  data_counts <- plotCounts(dds, gene = gene_id, intgroup = c("genotype", "infection"), returnData = TRUE)
  
  # Plot
  p_counts <- ggplot(data_counts, aes(x = interaction(genotype, infection), y = count, color = genotype)) +
    geom_point(position = position_jitter(w = 0.1, h = 0), size = 3) + 
    stat_summary(fun = mean, geom = "line", aes(group = 1), color = "gray", alpha = 0.5, linewidth = 1) +
    scale_y_log10() + 
    labs(title = paste("Expression of", gene_sym),
         subtitle = "Interaction Effect (Normalized Counts)",
         y = "Normalized Counts (log10)",
         x = "Condition") +
    theme_bw(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
      axis.text.y = element_text(size = 20),                           
      axis.title  = element_text(size = 24),                           
      plot.title  = element_text(size = 26, face = "bold"),            
      plot.subtitle = element_text(size = 22))
  
  print(p_counts)
  
  # Save Plot
  ggsave(file.path(base_dir, paste0("DE_Counts_", gene_sym, ".png")), plot = p_counts, width = 10, height = 8)
}


# ==============================================================================
# 10. FUNCTIONAL ENRICHMENT (GO) - SPLIT ANALYSIS (Strict Filter)
# ==============================================================================
# Analyze Up and Down genes separately. 
#Criteria: padj < 0.05 AND |Log2FC| > 1.

cat("Starting Split GO Analysis (Strict: |Log2FC| > 1)...\n")

# Define universe (all measured genes) excluding NA genes
universe_genes <- rownames(res_interaction)[!is.na(res_interaction$pvalue)]


# ------------------------------------------------------------------------------
# A. ANALYSIS OF DOWN-REGULATED GENES (Strictly < -1)
# ------------------------------------------------------------------------------
# These correspond to the BLUE dots in the Volcano Plot.
# Logic: take 'sig_genes' (already padj<0.05) and keep only Log2FC < -1.

genes_down_strict <- rownames(sig_genes[sig_genes$log2FoldChange < -1.0, ])

cat("Analyzing", length(genes_down_strict), "Down-regulated genes (Blue)...\n")

ego_down <- enrichGO(gene          = genes_down_strict,
                     universe      = universe_genes,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",     
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)

if(!is.null(ego_down) && nrow(ego_down) > 0) {
  
  # Barplot (Blue Theme)
  p_bar_down <- barplot(ego_down, showCategory=10) + 
    ggtitle("Top Enriched Terms (Down-regulated, Log2FC < -1)") +
    theme(plot.title = element_text(color = "navy", face = "bold", size = 24),
          axis.title = element_text(size = 24),                                 
          axis.text  = element_text(size = 16),                                 
          legend.title = element_text("p.adjust", size = 16),
          legend.text = element_text(size = 13))
  
  ggsave(file.path(base_dir, "GO_Barplot_Down.png"), p_bar_down, width = 13, height = 9)
  write.csv(as.data.frame(ego_down), file.path(base_dir, "GO_Down.csv"))
  
  print(p_bar_down)

} else {
  message("No GO terms found for strict Down-regulated genes.")
}


# ------------------------------------------------------------------------------
# B. ANALYSIS OF UP-REGULATED GENES (Strictly > 1)
# ------------------------------------------------------------------------------
# These correspond to the RED dots in the Volcano Plot.
# Logic: take 'sig_genes' (already padj<0.05) and keep only Log2FC > 1.

genes_up_strict <- rownames(sig_genes[sig_genes$log2FoldChange > 1.0, ])

cat("Analyzing", length(genes_up_strict), "Up-regulated genes (Red)...\n")

ego_up <- enrichGO(gene          = genes_up_strict,
                   universe      = universe_genes,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)

if(!is.null(ego_up) && nrow(ego_up) > 0) {
  
  # Barplot (Red Theme)
  p_bar_up <- barplot(ego_up, showCategory=10) + 
    ggtitle("Top Enriched Terms (Up-regulated, Log2FC > 1)") +
    theme(plot.title = element_text(color = "darkred", face = "bold", size = 24),
          axis.title = element_text(size = 24),                               
          axis.text  = element_text(size = 33),                                
          legend.title = element_text("p.adjust", size = 16),
          legend.text = element_text(size = 13))
  
  ggsave(file.path(base_dir, "GO_Barplot_Up.png"), p_bar_up, width = 13, height = 9)
  write.csv(as.data.frame(ego_up), file.path(base_dir, "GO_Up.csv"))
  
  print(p_bar_up)

} else {
  message("No GO terms found for strict Up-regulated genes.")
}

cat("Strict GO Analysis finished.\n")



# ------------------------------------------------------------------------------
# 11. Export Results Tables
# ------------------------------------------------------------------------------
write.csv(as.data.frame(res_interaction), file = file.path(base_dir, "DE_Interaction_Full_Results.csv"))
write.csv(as.data.frame(sig_genes), file = file.path(base_dir, "DE_Interaction_Significant_Genes.csv"))

message("Analysis complete. All outputs saved to: ", base_dir)