library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(GSVA)
library(GSVAdata)
library(SummarizedExperiment)
library(pheatmap)
library(RColorBrewer)
library(ashr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler) 
library(enrichR)
library(MASS)
library(VennDiagram)
library(msigdbr)
library(ggrepel)
library(fgsea)
library(svglite)
library(reshape2)
################################################################################
# Load the data
################################################################################
counts <- read.csv("~/The Francis Crick Dropbox/Daqi Deng/2_Data & Results/202502 RNA_seq/raw counts/corrected_counts_combat.csv", row.names=1, )
counts <- read.csv("./data/PDO_GE/combat_corrected_ctx_rna.csv", row.names=1 )

# subset protein_coding genes
pcg <- read.csv("./data/PDO_GE/protein_coding_genes.csv", sep="")
pcg <- pcg$symbol
genes <- intersect(pcg, rownames(counts))
counts <- counts[genes, ]
counts <- as.matrix(counts)   
rm(pcg, genes)

meta <- read.csv("./data/PDO_GE/meta_rna.csv",row.names = 1)
meta <- meta[meta$expts == "BGE", ]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
################################################################################
# for all the samples - related to Figure 6D
################################################################################
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation
# rld <- rlog(dds, blind = FALSE) # Regularised log transformation
plotPCA(vsd, intgroup=c("ge_status", "Line"))
vsd_counts <- assay(vsd)
# Calculate variance for each gene
gene_variances <- apply(vsd_counts, 1, var)
# Get indices of top 500 variable genes
top_genes <- order(gene_variances, decreasing = TRUE)[1:500]
# Subset expression matrix
vsd_counts <- vsd_counts[top_genes, ]

pca <- prcomp(t(vsd_counts))
sdev <- pca$sdev
variance_explained <- (sdev^2) / sum(sdev^2) * 100
# Extract variance explained for PC1 and PC2
pc1_variance <- variance_explained[1]
pc2_variance <- variance_explained[2]
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], BAP1_status = meta$BAP1_status, Background = meta$background)

# PCA plot with ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = BAP1_status, shape = Background)) +
  geom_point(size = 4) +
  scale_color_manual(
    values = c(
      "WT" = "#07519c",       # blue
      "KO" = "#cc4e26"  # red
    )
  ) +
  scale_shape_manual(values = c("Tumour" = 17, "Normal" = 16)) +  # 17 = triangle, 16 = circle
  labs(
    x = paste("PC1 (", round(pc1_variance, 1), "%)", sep = ""), 
    y = paste("PC2 (", round(pc2_variance, 1), "%)", sep = "")
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  )


sampleCor <- cor(assay(vsd), method="pearson")
# Convert correlation to a distance matrix (1 - correlation)
sampleDist <- as.dist(1 - sampleCor)
# Generate heatmap using correlation as distance metric
pheatmap(sampleCor, 
         main = "Sample Correlation Heatmap")
################################################################################
# K1040T1_GB1 BKO vs K1040T1_GB1 
################################################################################
meta <- meta[meta$line == "K1040T1_GB1", ]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation

# DESeq2 analysis
dds <- DESeq(dds)
contrast = c("ge_status", "bko", "wt")
res <- results(dds, contrast = contrast, alpha = 0.05)
res_lfcsk <- lfcShrink(dds, contrast = contrast, res=res, type = "ashr")
degs <- as.data.frame(res_lfcsk)
write.csv(degs, "./output/csv/K1040T1_GB1_deg.csv")
################################################################################
# Generate a MA plot using K1040T1_GB1 BKO vs K1040T1_GB1. related to Figure S8E
################################################################################
# Create a data frame for ggplot
plot_data <- data.frame(
  gene = rownames(res_lfcsk),
  baseMean = res_lfcsk$baseMean,
  log2FoldChange = res_lfcsk$log2FoldChange,
  padj = res_lfcsk$padj
)
gene_to_label <- plot_data[plot_data$gene %in% c("COL4A1", "COL4A2", "MYC", "CCND1", "ATF3", "JUNB", "CDC20", "CENPE", "CENPF", "KIF4A", "FN1", "CCL2", "MT1X", "SERPINE1","TGFB1", "MKI67"), ]
plot_data$direction <- "NS"
plot_data$direction[plot_data$log2FoldChange > 0 & plot_data$padj < 0.05] <- "Up"
plot_data$direction[plot_data$log2FoldChange < 0 & plot_data$padj < 0.05] <- "Down"

# MA plot with ggplot2
a<- ggplot(plot_data, aes(x = baseMean, y = log2FoldChange)) +
  geom_point(aes(color = direction), alpha = 1) +
  scale_color_manual(
    values = c("Up" = "#cc4e26", "Down" = "#07519c", "NS" = "grey"),
    labels = c("Up" = "Upregulated", "Down" = "Downregulated", "NS" = "Not significant"),
    name = "Direction"
  ) +
  scale_x_log10() +  # Log scale for mean expression
  labs(x = "baseMean", 
       y = "Log2FC") +
  ylim(-5, 5) +
  geom_hline(yintercept = 1, color = "red", linetype = "dashed") +
  geom_hline(yintercept = -1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 10, color = "red", linetype = "dashed") +
  geom_text_repel(
    data = gene_to_label,
    aes(label = gene),
    fontface = "bold.italic",
    box.padding = 0.5,  # Increase spacing
    point.padding = 0.3,
    segment.color = 'black',  # Linking lines
    segment.size = 0.5,  # Line thickness
    segment.linetype = "dashed",  # Dashed lines
    min.segment.length = 0,  # Ensures all labels have segments
    force = 2,  # Increase repulsion for better label placement
    direction = "both"  # Allow flexible movement
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # Add box/frame
    panel.grid = element_blank()  # Remove all gridlines
  )
ggsave(a, device = "svg", filename = "./output/figures/MA_plot_K1040_GE1.svg", width = 10, height = 6, dpi = 300)

################################################################################
# K1040T1_GB2 BKO vs K1040T1_GB2 
################################################################################
meta <- meta[meta$line == "K1040T1_GB2", ]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation

# DESeq2 analysis
dds <- DESeq(dds)
contrast = c("ge_status", "bko", "wt")
res <- results(dds, contrast = contrast, alpha = 0.05)
res_lfcsk <- lfcShrink(dds, contrast = contrast, res=res, type = "ashr")
degs <- as.data.frame(res_lfcsk)
write.csv(degs, "./output/csv/K1040T1_GB2_deg.csv")
################################################################################
# K1088 (CVBKO) vs K1088 (CVKO) 
################################################################################
meta <- meta[meta$background == "Normal",]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation

# DESeq2 analysis
dds <- DESeq(dds)
contrast = c("ge_status", "cvbko", "cvko")
res <- results(dds, contrast = contrast, alpha = 0.05)
res_lfcsk <- lfcShrink(dds, contrast = contrast, res=res, type = "ashr")
degs <- as.data.frame(res_lfcsk)
write.csv(degs, "./output/csv/K1088N1_GN1_deg_all.csv")
################################################################################
# PDOT (BKO) vs PDOT all
################################################################################
meta <- meta[meta$background == "Tumour", ]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation

# DESeq2 analysis
dds <- DESeq(dds)
contrast = c("ge_status", "bko", "wt")
res <- results(dds, contrast = contrast, alpha = 0.05)
res_lfcsk <- lfcShrink(dds, contrast = contrast, res=res, type = "ashr")
degs <- as.data.frame(res_lfcsk)
write.csv(degs, "./output/csv/K_PDOT_bko_deg.csv")
################################################################################
# Load all the DEG results
################################################################################
K1040_GB1_deg <- read.csv("./output/csv/K1040T1_GB1_deg.csv", row.names = 1)
K1040_GB2_deg <- read.csv("./output/csv/K1040T1_GB2_deg.csv", row.names = 1)
K1088_GN1_deg_all <- read.csv("./output/csv/K1088N1_GN1_deg_all.csv", row.names = 1)
K_PDOT_deg <- read.csv("./output/csv/K_PDOT_bko_deg.csv", row.names = 1)

################################################################################
# Generate an enrichment plot with K1088 GN DEGs and BAP1-specific genes. related to Figure S8C
################################################################################
BV <- read.csv("./data/PDO_DEG/BAP1vsN_all_markers.csv", row.names = 1)
PV <- read.csv("./data/PDO_DEG/PBRM1vsN_all_markers.csv", row.names = 1)
Buni <- setdiff(rownames(BV[BV$avg_log2FC > 0 & BV$p_val_adj < 0.01, ]), rownames(PV[PV$avg_log2FC > 0 & PV$p_val_adj < 0.01, ]))
Buni <- intersect(Buni, rownames(BV[BV$avg_log2FC>1, ]))
Buni <- list(BAP1_specific_genes = Buni)

# Create a named vector of logFC for ranking
ranked_genes <- K1088_GN1_deg_all$log2FoldChange
names(ranked_genes) <- rownames(K1088_GN1_deg_all)
# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Run fgsea
fgsea_results <- fgsea(
  pathways = Buni, 
  stats = ranked_genes, 
  minSize = 15, 
  maxSize = 500
)

# enrichemnt plot
nes <- round(fgsea_results$NES, 2)
fdr <- signif(fgsea_results$padj, 2)
es_data <- plotEnrichmentData(pathway = Buni[["BAP1_specific_genes"]], stats = ranked_genes)
# Extract components as data.frames
curve_df <- as.data.frame(es_data$curve)  # columns: rank, ES
ticks_df <- as.data.frame(es_data$ticks)  # columns: rank, stat
stats_df <- as.data.frame(es_data$stats)  # full ranked stat

# Rename for clarity
colnames(curve_df) <- c("x", "ES")
colnames(ticks_df) <- c("x", "stat")
colnames(stats_df) <- c("x", "rank_metric")

EP <- ggplot(curve_df, aes(x = x, y = ES)) +
  geom_line(color = "darkgreen", size = 1) + # Enrichment score curve
  geom_tile(
    data = stats_df,
    aes(x = x, y = -0.1, fill = rank_metric),
    height = 0.1,
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0, limits = c(-2, 2)) +  # Red–blue background OVERLAY (drawn later = on top)
  geom_linerange(
    data = ticks_df,
    aes(x = x, ymin = -0.15, ymax = -0.05),
    color = "black",
    alpha = 1,
    size = 0.1, 
    inherit.aes = FALSE
  ) +  # Tick marks for gene set members
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  ) +
  labs(
    y = "Enrichment Score (ES)",
    fill = "Log2FC",
    title = paste0("NES = ", nes, ", FDR = ", fdr))

################################################################################
# Generate an enrichment plot with PDOT & T-BKO and BAP1-specific genes. related to Figure S8D
################################################################################
# Create a named vector of logFC for ranking
ranked_genes <- K_PDOT_deg$log2FoldChange
names(ranked_genes) <- rownames(K_PDOT_deg)
# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Run fgsea
fgsea_results <- fgsea(
  pathways = Buni, 
  stats = ranked_genes, 
  minSize = 15, 
  maxSize = 500
)

# enrichemnt plot
nes <- round(fgsea_results$NES, 2)
fdr <- signif(fgsea_results$padj, 2)
es_data <- plotEnrichmentData(pathway = Buni[["BAP1_specific_genes"]], stats = ranked_genes)
# Extract components as data.frames
curve_df <- as.data.frame(es_data$curve)  # columns: rank, ES
ticks_df <- as.data.frame(es_data$ticks)  # columns: rank, stat
stats_df <- as.data.frame(es_data$stats)  # full ranked stat

# Rename for clarity
colnames(curve_df) <- c("x", "ES")
colnames(ticks_df) <- c("x", "stat")
colnames(stats_df) <- c("x", "rank_metric")

EP <- ggplot(curve_df, aes(x = x, y = ES)) +
  geom_line(color = "darkgreen", size = 1) + # Enrichment score curve
  geom_tile(
    data = stats_df,
    aes(x = x, y = -0.1, fill = rank_metric),
    height = 0.1,
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0, limits = c(-2, 2)) +  # Red–blue background OVERLAY (drawn later = on top)
  geom_linerange(
    data = ticks_df,
    aes(x = x, ymin = -0.15, ymax = -0.05),
    color = "black",
    alpha = 1,
    size = 0.1, 
    inherit.aes = FALSE
  ) +  # Tick marks for gene set members
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  ) +
  labs(
    y = "Enrichment Score (ES)",
    fill = "Log2FC",
    title = paste0("NES = ", nes, ", FDR = ", fdr))

################################################################################
# Gene set enrichemnt analysis
################################################################################
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")

hallmark_pathways <- hallmark_genesets %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Create a named vector of logFC for ranking
deg <- K_PDOT_deg
ranked_genes <- deg$log2FoldChange
names(ranked_genes) <- rownames(deg)
# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Run fgsea
fgsea_results <- fgsea(
  pathways = hallmark_pathways, 
  stats = ranked_genes, 
  minSize = 15, 
  maxSize = 500
)

################################################################################
# Generate an enrichment plot with PDOT DEGs and hallmark cell cycle 
################################################################################
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_pathways <- hallmark_genesets %>%
  split(x = .$gene_symbol, f = .$gs_name)
E2F_target <- list(HALLMARK_E2F_TARGETS = hallmark_pathways[["HALLMARK_E2F_TARGETS"]])

# Create a named vector of logFC for ranking
ranked_genes <- K_PDOT_deg$log2FoldChange
names(ranked_genes) <- rownames(K_PDOT_deg)
# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Run fgsea
fgsea_results <- fgsea(
  pathways = E2F_target, 
  stats = ranked_genes, 
  minSize = 15, 
  maxSize = 500
)

# enrichemnt plot
nes <- round(fgsea_results$NES, 2)
fdr <- signif(fgsea_results$padj, 2)
es_data <- plotEnrichmentData(pathway = E2F_target[["HALLMARK_E2F_TARGETS"]], stats = ranked_genes)
# Extract components as data.frames
curve_df <- as.data.frame(es_data$curve)  # columns: rank, ES
ticks_df <- as.data.frame(es_data$ticks)  # columns: rank, stat
stats_df <- as.data.frame(es_data$stats)  # full ranked stat

# Rename for clarity
colnames(curve_df) <- c("x", "ES")
colnames(ticks_df) <- c("x", "stat")
colnames(stats_df) <- c("x", "rank_metric")

EP <- ggplot(curve_df, aes(x = x, y = ES)) +
  geom_line(color = "darkgreen", size = 1) + # Enrichment score curve
  geom_tile(
    data = stats_df,
    aes(x = x, y = -0.8, fill = rank_metric),
    height = 0.1,
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0, limits = c(-2, 2)) +  # Red–blue background OVERLAY (drawn later = on top)
  geom_linerange(
    data = ticks_df,
    aes(x = x, ymin = -0.85, ymax = -0.75),
    color = "black",
    alpha = 1,
    size = 0.1, 
    inherit.aes = FALSE
  ) +  # Tick marks for gene set members
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  ) +
  labs(
    y = "Enrichment Score (ES)",
    fill = "Log2FC",
    title = paste0("NES = ", nes, ", FDR = ", fdr))

################################################################################
# Generate an enrichment plot with PDOT DEGs and hallmark TGFB 
################################################################################
# get the top BAP1 genes FC > 2
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_pathways <- hallmark_genesets %>%
  split(x = .$gene_symbol, f = .$gs_name)
TGFB <- list(HALLMARK_TGF_BETA_SIGNALING = hallmark_pathways[["HALLMARK_TGF_BETA_SIGNALING"]])

# Create a named vector of logFC for ranking
ranked_genes <- K_PDOT_deg$log2FoldChange
names(ranked_genes) <- rownames(K_PDOT_deg)
# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Run fgsea
fgsea_results <- fgsea(
  pathways = TGFB, 
  stats = ranked_genes, 
  minSize = 15, 
  maxSize = 500
)

# enrichemnt plot
nes <- round(fgsea_results$NES, 2)
fdr <- signif(fgsea_results$padj, 2)
es_data <- plotEnrichmentData(pathway = TGFB[["HALLMARK_TGF_BETA_SIGNALING"]], stats = ranked_genes)
# Extract components as data.frames
curve_df <- as.data.frame(es_data$curve)  # columns: rank, ES
ticks_df <- as.data.frame(es_data$ticks)  # columns: rank, stat
stats_df <- as.data.frame(es_data$stats)  # full ranked stat

# Rename for clarity
colnames(curve_df) <- c("x", "ES")
colnames(ticks_df) <- c("x", "stat")
colnames(stats_df) <- c("x", "rank_metric")

EP <- ggplot(curve_df, aes(x = x, y = ES)) +
  geom_line(color = "darkgreen", size = 1) + # Enrichment score curve
  geom_tile(
    data = stats_df,
    aes(x = x, y = -0.1, fill = rank_metric),
    height = 0.1,
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0, limits = c(-2, 2)) +  # Red–blue background OVERLAY (drawn later = on top)
  geom_linerange(
    data = ticks_df,
    aes(x = x, ymin = -0.15, ymax = -0.05),
    color = "black",
    alpha = 1,
    size = 0.1, 
    inherit.aes = FALSE
  ) +  # Tick marks for gene set members
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  ) +
  labs(
    y = "Enrichment Score (ES)",
    fill = "Log2FC",
    title = paste0("NES = ", nes, ", FDR = ", fdr))

################################################################################
# Generate an enrichment plot with PDOT DEGs and hallmark IL2_STAT5_SIGNALING
################################################################################
# get the top BAP1 genes FC > 2
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_pathways <- hallmark_genesets %>%
  split(x = .$gene_symbol, f = .$gs_name)
IL2_STAT5_SIGNALING <- list(HALLMARK_IL2_STAT5_SIGNALING = hallmark_pathways[["HALLMARK_IL2_STAT5_SIGNALING"]])

# Create a named vector of logFC for ranking
ranked_genes <- K_PDOT_deg$log2FoldChange
names(ranked_genes) <- rownames(K_PDOT_deg)
# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Run fgsea
fgsea_results <- fgsea(
  pathways = IL2_STAT5_SIGNALING, 
  stats = ranked_genes, 
  minSize = 15, 
  maxSize = 500
)

# enrichemnt plot
nes <- round(fgsea_results$NES, 2)
fdr <- signif(fgsea_results$padj, 2)
es_data <- plotEnrichmentData(pathway = IL2_STAT5_SIGNALING[["HALLMARK_IL2_STAT5_SIGNALING"]], stats = ranked_genes)
# Extract components as data.frames
curve_df <- as.data.frame(es_data$curve)  # columns: rank, ES
ticks_df <- as.data.frame(es_data$ticks)  # columns: rank, stat
stats_df <- as.data.frame(es_data$stats)  # full ranked stat

# Rename for clarity
colnames(curve_df) <- c("x", "ES")
colnames(ticks_df) <- c("x", "stat")
colnames(stats_df) <- c("x", "rank_metric")

EP <- ggplot(curve_df, aes(x = x, y = ES)) +
  geom_line(color = "darkgreen", size = 1) + # Enrichment score curve
  geom_tile(
    data = stats_df,
    aes(x = x, y = -0.2, fill = rank_metric),
    height = 0.1,
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0, limits = c(-2, 2)) +  # Red–blue background OVERLAY (drawn later = on top)
  geom_linerange(
    data = ticks_df,
    aes(x = x, ymin = -0.25, ymax = -0.15),
    color = "black",
    alpha = 1,
    size = 0.1, 
    inherit.aes = FALSE
  ) +  # Tick marks for gene set members
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  ) +
  labs(
    y = "Enrichment Score (ES)",
    fill = "Log2FC",
    title = paste0("NES = ", nes, ", FDR = ", fdr))

################################################################################
# Generate an enrichment plot with PDOT DEGs and GO Chromosome segregation 
################################################################################
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
hallmark_pathways <- hallmark_genesets %>%
  split(x = .$gene_symbol, f = .$gs_name)
CIN <- list(GOBP_CHROMOSOME_SEGREGATION = hallmark_pathways[["GOBP_CHROMOSOME_SEGREGATION"]])

# Create a named vector of logFC for ranking
ranked_genes <- K_PDOT_deg$log2FoldChange
names(ranked_genes) <- rownames(K_PDOT_deg)
# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Run fgsea
fgsea_results <- fgsea(
  pathways = CIN, 
  stats = ranked_genes, 
  minSize = 15, 
  maxSize = 500
)

# enrichemnt plot
nes <- round(fgsea_results$NES, 2)
fdr <- signif(fgsea_results$padj, 2)
es_data <- plotEnrichmentData(pathway = CIN[["GOBP_CHROMOSOME_SEGREGATION"]], stats = ranked_genes)
# Extract components as data.frames
curve_df <- as.data.frame(es_data$curve)  # columns: rank, ES
ticks_df <- as.data.frame(es_data$ticks)  # columns: rank, stat
stats_df <- as.data.frame(es_data$stats)  # full ranked stat

# Rename for clarity
colnames(curve_df) <- c("x", "ES")
colnames(ticks_df) <- c("x", "stat")
colnames(stats_df) <- c("x", "rank_metric")

EP <- ggplot(curve_df, aes(x = x, y = ES)) +
  geom_line(color = "darkgreen", size = 1) + # Enrichment score curve
  geom_tile(
    data = stats_df,
    aes(x = x, y = -0.7, fill = rank_metric),
    height = 0.1,
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkred", midpoint = 0, limits = c(-2, 2)) +  # Red–blue background OVERLAY (drawn later = on top)
  geom_linerange(
    data = ticks_df,
    aes(x = x, ymin = -0.75, ymax = -0.65),
    color = "black",
    alpha = 1,
    size = 0.1, 
    inherit.aes = FALSE
  ) +  # Tick marks for gene set members
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  ) +
  labs(
    y = "Enrichment Score (ES)",
    fill = "Log2FC",
    title = paste0("NES = ", nes, ", FDR = ", fdr))


