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
setwd("~/Downloads/PBRM1_BAP1_KidneyCancer")
################################################################################
# Load the data - this needs to run before subsetting for each line
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
meta <- meta[meta$expts == "PGE", ]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
################################################################################
# for all the samples - related to Figure 5B
################################################################################
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation
# rld <- rlog(dds, blind = FALSE) # Regularised log transformation
plotPCA(vsd, intgroup=c("ge_status", "batch"))
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
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PBRM1_Status = meta$PBRM1_status, Background = meta$background)

# PCA plot with ggplot2
ggplot(pca_data, aes(x = PC1, y = PC2, color = PBRM1_Status, shape = Background)) +
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
# K1030T1_GB1 PKO vs K1030T1_GB1
################################################################################
meta <- meta[meta$line == "K1030T1_GB1", ]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation

# DESeq2 analysis
dds <- DESeq(dds)
contrast = c("ge_status", "pko", "wt")
res <- results(dds, contrast = contrast, alpha = 0.05)
res_lfcsk <- lfcShrink(dds, contrast = contrast, res=res, type = "ashr")
degs <- as.data.frame(res_lfcsk)
write.csv(degs, "./output/csv/K1030T1_GB1_deg_all.csv")
################################################################################
# K1030T2_GB1 PKO vs K1030T2_GB1
################################################################################
meta <- meta[meta$line == "K1030T2_GB1", ]
meta <- meta[meta$batch == 2, ]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation

# DESeq2 analysis
dds <- DESeq(dds)
contrast = c("ge_status", "pko", "wt")
res <- results(dds, contrast = contrast, alpha = 0.05)
res_lfcsk <- lfcShrink(dds, contrast = contrast, res=res, type = "ashr")
degs <- as.data.frame(res_lfcsk)
write.csv(degs, "./output/csv/K1030T2_GB1_deg_all.csv")
################################################################################
# K1063T1_GB1 PKO vs K1063T1_GB1 
################################################################################
meta <- meta[meta$line == "K1063T1_GB1", ]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation

# DESeq2 analysis
dds <- DESeq(dds)
contrast = c("ge_status", "pko", "wt")
res <- results(dds, contrast = contrast, alpha = 0.05)
res_lfcsk <- lfcShrink(dds, contrast = contrast, res=res, type = "ashr")
degs <- as.data.frame(res_lfcsk)
write.csv(degs, "./output/csv/K1063T1_GB1_deg_all.csv")
################################################################################
# PDON (PKO) vs PDON 
################################################################################
meta <- meta[meta$background == "Normal", ]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation

# DESeq2 analysis
dds <- DESeq(dds)
contrast = c("ge_status", "pko", "wt")
res <- results(dds, contrast = contrast, alpha = 0.05)
res_lfcsk <- lfcShrink(dds, contrast = contrast, res=res, type = "ashr")
degs <- as.data.frame(res_lfcsk)
write.csv(degs, "./output/csv/K_PDON_pko_deg_all.csv")
################################################################################
# PDOT (PKO) vs PDOT all
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
contrast = c("ge_status", "pko", "wt")
res <- results(dds, contrast = contrast, alpha = 0.05)
res_lfcsk <- lfcShrink(dds, contrast = contrast, res=res, type = "ashr")
degs <- as.data.frame(res_lfcsk)
write.csv(degs, "./output/csv/K_PDOT_pko_deg_all.csv")
################################################################################
# Load all the DEG results for analysis
################################################################################
K_PDON_deg <- read.csv("./output/csv/K_PDON_pko_deg_all.csv", row.names = 1)
K_PDOT_deg <- read.csv("./output/csv/K_PDOT_pko_deg_all.csv", row.names = 1)
K1030_GB1_deg <- read.csv("./output/csv/K1030T1_GB1_deg_all.csv", row.names = 1)
K1030_GB2_deg <- read.csv("./output/csv/K1030T2_GB1_deg_all.csv", row.names = 1)
K1063_GB1_deg <- read.csv("./output/csv/K1063T1_GB1_deg_all.csv", row.names = 1)

################################################################################
# Generate an enrichment plot with K_PDON_deg  DEGs and BAP1-specific genes related to Figure S7E
################################################################################
# get the top BAP1 genes FC > 2
BV <- read.csv("./data/PDO_DEG/BAP1vsN_all_markers.csv", row.names = 1)
PV <- read.csv("./data/PDO_DEG/PBRM1vsN_all_markers.csv", row.names = 1)
Puni <- setdiff(rownames(PV[PV$avg_log2FC > 0 & PV$p_val_adj < 0.01, ]), rownames(BV[BV$avg_log2FC > 0 & BV$p_val_adj < 0.01, ]))
Puni <- list(PBRM1_specific_genes = Puni)

# Create a named vector of logFC for ranking
ranked_genes <- K_PDON_deg$log2FoldChange
names(ranked_genes) <- rownames(K_PDON_deg)
# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Run fgsea
fgsea_results <- fgsea(
  pathways = Puni, 
  stats = ranked_genes, 
  minSize = 15, 
  maxSize = 1000
)

# enrichemnt plot
nes <- round(fgsea_results$NES, 2)
fdr <- signif(fgsea_results$padj, 2)
es_data <- plotEnrichmentData(pathway = Puni[["PBRM1_specific_genes"]], stats = ranked_genes)
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
ggsave("./results/K_PDON_deg_PBRM1_specific_gene_enrichment_plot.svg", 
       plot = EP, 
       width = 4, 
       height = 2.5, 
       device = "svg")
################################################################################
# Generate an enrichment plot with K_PDOT_deg  DEGs and BAP1-specific genes related to Figure S7F
################################################################################

# Create a named vector of logFC for ranking
ranked_genes <- K_PDOT_deg$log2FoldChange
names(ranked_genes) <- rownames(K_PDOT_deg)
# Sort in decreasing order
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Run fgsea
fgsea_results <- fgsea(
  pathways = Puni, 
  stats = ranked_genes, 
  minSize = 15, 
  maxSize = 1000
)

# enrichemnt plot
nes <- round(fgsea_results$NES, 2)
fdr <- signif(fgsea_results$padj, 2)
es_data <- plotEnrichmentData(pathway = Puni[["PBRM1_specific_genes"]], stats = ranked_genes)
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
ggsave("./results//K_PDOT_deg_PBRM1_specific_gene_enrichment_plot.svg", 
       plot = EP, 
       width = 4, 
       height = 2.5, 
       device = "svg")
################################################################################
# PRODUCE A RANK PLOT USING K1030_GB2_deg. related to Figure 5C
################################################################################
K1030_GB2_deg <- K1030_GB2_deg[order(K1030_GB2_deg$log2FoldChange, decreasing = TRUE), ]
K1030_GB2_deg$rank <- 1:nrow(K1030_GB2_deg)
ggplot(K1030_GB2_deg, aes(x = rank, y = log2FoldChange)) +
  geom_line(color = "grey70") +
  geom_point(size = 0.7, alpha = 0.8) +
  labs(title = "Ranked Differential Expression Plot",
       x = "Gene Rank",
       y = "log2 Fold Change") +
  theme_minimal()
genes_to_annotate <- c("WNT7B", "WNT5B", "IGF2", "FGF9", "IGF1")

K1030_GB2_deg$label <- ifelse(rownames(K1030_GB2_deg) %in% genes_to_annotate, rownames(K1030_GB2_deg), NA)

ggplot(K1030_GB2_deg, aes(x = rank, y = log2FoldChange)) +
  geom_line(color = "grey80") +
  geom_point(size = 0.7, alpha = 0.8) +
  geom_text_repel(
    aes(label = label),
    size = 6,
    max.overlaps = 10,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "black",       # line color
    segment.size = 1,            # line thickness
    min.segment.length = 0         # always draw line
  ) +
  labs(x = "Gene Rank",
       y = "log2 Fold Change") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  )
################################################################################
# PRODUCE A VENN DIAGRAM. related to Figure 5F
################################################################################

genes_A <- rownames(K1030_GB1_deg[K1030_GB1_deg$log2FoldChange > 0, ])
genes_B <- rownames(K1030_GB2_deg[K1030_GB2_deg$log2FoldChange > 0, ])
genes_C <- rownames(K1063_GB1_deg[K1063_GB1_deg$log2FoldChange > 0, ])

# Draw 3-set Venn
venn.plot <- venn.diagram(
  x = list(
    K1030_GB1 = genes_A,
    K1030_GB2 = genes_B,
    K1063_GB1 = genes_C
  ),
  filename = NULL,        # Don't save to file
  fill = c("red", "blue", "green"),
  alpha = 0.4,
  cex = 1.5,
  cat.cex = 1.3,
  cat.pos = 0,
  margin = 0.1
)
grid::grid.draw(venn.plot)
################################################################################
# scatter plot of PDON and PDOT related to Figure 5G
################################################################################
K_PDON_deg$gene <- rownames(K_PDON_deg)
K_PDOT_deg$gene <- rownames(K_PDOT_deg)

merged_deg <- merge(K_PDON_deg, K_PDOT_deg, by = "gene", suffixes = c("_A", "_B"))

merged_deg$quadrant <- with(merged_deg, ifelse(
  log2FoldChange_A > 0 & log2FoldChange_B > 0, "Up-Up",
  ifelse(log2FoldChange_A < 0 & log2FoldChange_B < 0, "Down-Down",
         ifelse(log2FoldChange_A > 0 & log2FoldChange_B < 0, "Up-Down", "Down-Up"))
))

ggplot(merged_deg, aes(x = log2FoldChange_A, y = log2FoldChange_B, color = quadrant)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  scale_x_continuous(
    limits = c(-5,5),
    breaks = c(-5, 0, 5)
  ) +
  scale_y_continuous(
    limits = c(-4,4),
    breaks = c(-4, 0, 4)
  ) +
  coord_fixed() +
  scale_color_manual(values = c("Up-Up" = "#cc4e26", "Down-Down" = "#07519c", 
                                "Up-Down" = "grey", "Down-Up" = "grey")) +
  theme_minimal() +
  theme(legend.position = "none") +  # <- This removes the legend
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
  ) +
  labs(
    x = "Log2FC in Normal Background",   # <- customised x-axis label
    y = "Log2FC in Tumour Backgroud"    # <- customised y-axis label
  )
################################################################################
# GENERATE A z-score heatmap for the constant genes. related to Figure 5H
################################################################################
const_list <- merged_deg[merged_deg$quadrant %in% c("Up-Up", "Down-Down"), ]
# meta <- meta[meta$batch == , ]
vst_mat <- as.data.frame(assay(vsd))
vst_mat <- vst_mat[const_list$gene, ]
meta <- meta[!(rownames(meta) %in% c("K1030_GB210", "K1030_GB213", "K1030_GB223_TRI")), ]
vst_mat <- vst_mat[,rownames(meta)]
meta$group <- paste0(meta$background, "_", meta$ge_status)
meta$sample <- rownames(meta)
group_levels <- unique(meta$group)

# Result: genes x groups matrix
vst_group_means <- sapply(group_levels, function(grp) {
  samples_in_grp <- meta$sample[meta$group == grp]
  rowMeans(vst_mat[, samples_in_grp, drop = FALSE])
})
vst_z <- t(scale(t(vst_group_means)))  # genes x groups
vst_z <- vst_z[, c("Normal_wt", "Normal_pko", "Tumour_wt", "Tumour_pko")]

vst_z_df <- melt(vst_z)
colnames(vst_z_df) <- c("Gene", "Group", "Z-score")

genes_to_label <- c("JAK2", "STAT4", "SYNE2", "TGFBI", "TSC1", "YBX1", "ANGPT2", "ARF1")  # example gene names
row_labels <- ifelse(rownames(vst_z) %in% genes_to_label,
                     rownames(vst_z), "")

pheatmap(
  vst_z,
  cluster_rows = TRUE,       # cluster genes
  cluster_cols = F,      # keep group order fixed
  treeheight_row = 0,      # hide row dendrogram
  scale = "none",            # already Z-scored
  show_rownames = T,      # set FALSE if too many genes
  show_colnames = TRUE,
  labels_row = row_labels,  # custom labels
  color = colorRampPalette(c("#07519c", "white", "#cc4e26"))(100),
  fontsize_row = 6,          # adjust for clarity
  border_color = NA,
)
################################################################################
# Create a MA plot related to Figure S7G
################################################################################
meta <- meta[meta$line == "K1030T2_GB1", ]
meta <- meta[meta$batch == 2, ]
counts <- counts[,rownames(meta)]
all(rownames(meta) == colnames(counts)) # Should return TRUE
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData = meta,
                              design = ~ ge_status)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- vst(dds, blind = FALSE)  # Variance-stabilizing transformation

# DESeq2 analysis
dds <- DESeq(dds)
contrast = c("ge_status", "pko", "wt")
res <- results(dds, contrast = contrast, alpha = 0.05)
res_lfcsk <- lfcShrink(dds, contrast = contrast, res=res, type = "ashr")
plot_data <- data.frame(
  gene = rownames(res_lfcsk),
  baseMean = res_lfcsk$baseMean,
  log2FoldChange = res_lfcsk$log2FoldChange,
  padj = res_lfcsk$padj
)
gene_to_label <- plot_data[plot_data$gene %in% c("IGF1", "LOX", "IGF2", "IGFBP7", "POU5F1", "MMP2", "WNT5B", "WNT7B", "ANPEP", "BMP1", "ELK3", "FGF9", "FN1", "IGFBP6", "ITGAV", "ITGB2", "ITGB3", "KLF4", "CD24", "PROM1", "IGFBP2", "CCL2", "ANGPT2", "CDH1"), ]
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
ggsave(a, device = "svg", filename = "./output/figures/MA_plot_K1030_GE5.svg", width = 10, height = 6, dpi = 300)
################################################################################

