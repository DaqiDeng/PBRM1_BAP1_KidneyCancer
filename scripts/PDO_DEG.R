library(dplyr)
library(tidyr)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggsankey)
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
################################################################################
################################################################################
# Load the data
################################################################################
setwd("~/Downloads/PBRM1_BAP1_KidneyCancer/")

file_list <- list.files(path = "./data/PDO_DEG", pattern = "*.csv", full.names = T)
for (file in file_list) {
  df_name <- gsub(".*/|\\.csv$", "", file)
  assign(df_name, read.csv(file, row.names = 1))
}
remove(file_list, file, df_name)
background <- background$x

################################################################################
# to define and get different gene sets - venn diagram for Figure 2E&F
################################################################################
BupN <- rownames(BAP1vsN_all_markers[BAP1vsN_all_markers$avg_log2FC > 0, ])
BupV <- rownames(BAP1vsV_all_markers[BAP1vsV_all_markers$avg_log2FC > 0, ])
PupN <- rownames(PBRM1vsN_all_markers[PBRM1vsN_all_markers$avg_log2FC > 0, ])
PupV <- rownames(PBRM1vsV_all_markers[PBRM1vsV_all_markers$avg_log2FC > 0, ])

Bunique_baseN <- setdiff(BupN, PupN)
Punique_baseN <- setdiff(PupN, BupN)
shared_baseN <- intersect(BupN, PupN)

Bunique_baseV <- setdiff(BupV, PupV)
Punique_baseV <- setdiff(PupV, BupV)
shared_baseV <- intersect(BupV, PupV)

BaseN <- venn.diagram(
  x = list(Set1 = BupN, Set2 = PupN), 
  filename = NULL,
  category.names = c("VB vs N\nDEGs", "VP vs N\nDEGs"),
  col = c("#F37E78", "#4189CB"),
  fill = c("#F37E78", "#4189CB"),
  alpha = 0.2,
  cex = 2,
  cat.cex = 2,  # Font size for category labels
  cat.pos = c(-20, 20),  # Adjust label positions
  cat.dist = 0.05  # Adjust label distance from the circle
)
grid.draw(BaseN)

BaseV <- venn.diagram(
  x = list(Set1 = BupV, Set2 = PupV), 
  filename = NULL,
  category.names = c("VB vs V\nDEGs", "VP vs V\nDEGs"),
  col = c("#F37E78", "#4189CB"),
  fill = c("#F37E78", "#4189CB"),
  alpha = 0.2,
  cex = 2,
  cat.cex = 2,  # Font size for category labels
  cat.pos = c(-20, 20),  # Adjust label positions
  cat.dist = 0.05  # Adjust label distance from the circle
)
grid.draw(BaseV)

################################################################################
# NES calulation - for Figure 2G 
################################################################################
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_pathways <- hallmark_genesets %>%
  split(x = .$gene_symbol, f = .$gs_name)

BN <- setNames(BAP1vsN_all_markers$avg_log2FC, rownames(BAP1vsN_all_markers))
BN <- sort(BN, decreasing = TRUE)
BN_fgsea <- fgsea(
  pathways = hallmark_pathways,  # Your hallmark gene sets
  stats = BN,           # Ranked gene list
  minSize = 15,                  # Minimum size of a gene set to be tested (default typically between 15-25)
  maxSize = 500,                 # Maximum size of a gene set to be tested
)

PN <- setNames(PBRM1vsN_all_markers$avg_log2FC, rownames(PBRM1vsN_all_markers))
PN <- sort(PN, decreasing = TRUE)
PN_fgsea <- fgsea(
  pathways = hallmark_pathways,  # Your hallmark gene sets
  stats = PN,           # Ranked gene list
  minSize = 15,                  # Minimum size of a gene set to be tested (default typically between 15-25)
  maxSize = 500,                 # Maximum size of a gene set to be tested
)

BV <- setNames(BAP1vsV_all_markers$avg_log2FC, rownames(BAP1vsV_all_markers))
BV <- sort(BV, decreasing = TRUE)
BV_fgsea <- fgsea(
  pathways = hallmark_pathways,  # Your hallmark gene sets
  stats = BV,           # Ranked gene list
  minSize = 15,                  # Minimum size of a gene set to be tested (default typically between 15-25)
  maxSize = 500,                 # Maximum size of a gene set to be tested
)

PV <- setNames(PBRM1vsV_all_markers$avg_log2FC, rownames(PBRM1vsV_all_markers))
PV <- sort(PV, decreasing = TRUE)
PV_fgsea <- fgsea(
  pathways = hallmark_pathways,  # Your hallmark gene sets
  stats = PV,           # Ranked gene list
  minSize = 15,                  # Minimum size of a gene set to be tested (default typically between 15-25)
  maxSize = 500,                 # Maximum size of a gene set to be tested
)

fgsea_results_all <- bind_rows(
  BN_fgsea %>% mutate(gene_set = "BN"),
  PN_fgsea %>% mutate(gene_set = "PN"),
  BV_fgsea %>% mutate(gene_set = "BV"),
  PV_fgsea %>% mutate(gene_set = "PV")
)

fgsea_results_all <- fgsea_results_all %>%
  group_by(pathway) %>%  # Group by pathway name
  filter(any(padj <= 0.1)) %>%  # Keep only pathways with at least one padj <= 0.1
  ungroup()  # Ungroup for further analysis
fgsea_results_all$pathway <- gsub("HALLMARK_", "", fgsea_results_all$pathway)

# Bubble plot with customized axes and bubble aesthetics
# Bubble plot with customized axes, size, and color for NES
a <- ggplot(fgsea_results_all, aes(x = factor(gene_set, levels = c("BN", "PN", "BV", "PV")), y = reorder(pathway, NES), 
                              size = -log10(padj), color = NES)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +  # Color gradient for NES
  theme_minimal() +
  labs(title = "Bubble Plot of Gene Set Enrichment by Pathway",
       x = "Gene Set",
       y = "Pathway",
       color = "NES",
       size = "-log10padj") +
  theme(axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

ggsave(filename = "./output/figures/PDO_DEG_NES_bubble_plot.svg", plot = a, device = "svg", width = 8, height = 8)


################################################################################
# Volcano plots
################################################################################
log2FC_threshold <- 1
pvalue_threshold <- 0.01
BAP1vsPBRM1_markers$p_val_adj <- ifelse(BAP1vsPBRM1_markers$p_val_adj == 0, 1e-300, BAP1vsPBRM1_markers$p_val_adj)
BAP1vsPBRM1_markers$minusLog10Pvalue <- -log10(BAP1vsPBRM1_markers$p_val_adj)
BAP1vsPBRM1_markers$X <- rownames(BAP1vsPBRM1_markers)
gene_to_label <- BAP1vsPBRM1_markers[BAP1vsPBRM1_markers$X %in% c("VEGFA", "GRM8", "HAVCR1", "ANGPTL4", "GRIN2A", "HIF1A", "SOX6", "COL8A1", "TGFBI", "SERPINE1", "COL6A2", "PBX3", "FN1", "EPO", "GEM", "VEGFB", "COL4A5", "TMEM45A", "EPAS1", "COL6A1"),]


volcano_plot <- ggplot(BAP1vsPBRM1_markers, aes(x = avg_log2FC, y = minusLog10Pvalue)) +
  geom_point(aes(color = ifelse(p_val_adj > pvalue_threshold | abs(avg_log2FC) < log2FC_threshold, 'Non-significant', ifelse(p_val_adj < pvalue_threshold & avg_log2FC > log2FC_threshold & pct.1 < 0.2, 'B_pct<0.2', ifelse(p_val_adj < pvalue_threshold & avg_log2FC < -log2FC_threshold & pct.2 < 0.2, 'P_pct<0.2', ifelse(p_val_adj < pvalue_threshold & avg_log2FC > log2FC_threshold & pct.1 > 0.2, "BAP1_up", 'PBRM1_up')))), alpha = 0.7)) +
  scale_color_manual(values = c('Non-significant' = "grey90", 'B_pct<0.2' = "grey", 'P_pct<0.2' = "grey", 'BAP1_up' = "red", "PBRM1_up" = "red")) +
  geom_text_repel(data = gene_to_label, aes(label=X), 
                  box.padding = 0.35, point.padding = 0.5,
                  segment.color = 'grey50') +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Significant") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
volcano_plot
