library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggsankey)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(GSVA)
library(GSVAdata)
library(AnnotationDbi)
library(enrichR)
library(fgsea)
library(svglite)
setwd("~/Downloads/PBRM1_BAP1_KidneyCancer/")
################################################################################
# # Read in files
################################################################################
file_list <- list.files(path = "./data/scRNA/Markers/", pattern = "*.csv", full.names = T)
for (file in file_list) {
  df_name <- gsub(".*/|\\.csv$", "", file)
  assign(df_name, read.csv(file, row.names = 1))
}
remove(file_list, file, df_name)

################################################################################
# # Enrichemnt
################################################################################
# filter by percent expressed > 10%
BAP1_all_markers <- BAP1_all_std_markers[(BAP1_all_std_markers$avg_diff > 0 & BAP1_all_std_markers$pct.1 > 0.1) | (BAP1_all_std_markers$avg_diff < 0 & BAP1_all_std_markers$pct.2 > 0.1), ]
PBRM1_all_markers <- PBRM1_all_std_markers[(PBRM1_all_std_markers$avg_diff > 0 & PBRM1_all_std_markers$pct.1 > 0.1) | (PBRM1_all_std_markers$avg_diff < 0 & PBRM1_all_std_markers$pct.2 > 0.1), ]

BAP1_list <- rownames(BAP1_all_markers[BAP1_all_markers$avg_diff >0,])
PBRM1_list <- rownames(PBRM1_all_markers[PBRM1_all_markers$avg_diff >0,])

# Run MSigDB geneset enrichment
# retrieve all human gene sets - H: hallmark genesets
m_df <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_df, 2) %>% as.data.frame

em_BV_B <- enricher(BAP1_list, TERM2GENE=m_df, universe = background$x)

em_PV_P <- enricher(PBRM1_list, TERM2GENE=m_df, universe = background$x)

em_BV_B_short <- em_BV_B@result
em_PV_P_short <- em_PV_P@result

em_BV_B_short <- em_BV_B_short[em_BV_B_short$ID %in% c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_HYPOXIA"), ]

em_PV_P_short <- em_PV_P_short[em_PV_P_short$ID %in% c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_HYPOXIA"), ]

em_BV_B_short$Description <- substr(em_BV_B_short$Description, 10, 50)
em_PV_P_short$Description <- substr(em_PV_P_short$Description, 10, 50)

em_BV_B_short <- em_BV_B_short %>%
  mutate(
    log10padj = -log10(p.adjust),
    sig_group = ifelse(p.adjust < 0.05, "Significant", "Non-significant"),
    Description = factor(Description, c("EPITHELIAL_MESENCHYMAL_TRANSITION", "TNFA_SIGNALING_VIA_NFKB", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", "INFLAMMATORY_RESPONSE", "ANDROGEN_RESPONSE", "ESTROGEN_RESPONSE_EARLY", "HYPOXIA")),
    GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))
em_PV_P_short <- em_PV_P_short %>%
  mutate(
    log10padj = -log10(p.adjust),
    sig_group = ifelse(p.adjust < 0.05, "Significant", "Non-significant"),
    Description = factor(Description, c("EPITHELIAL_MESENCHYMAL_TRANSITION", "TNFA_SIGNALING_VIA_NFKB", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", "INFLAMMATORY_RESPONSE", "ANDROGEN_RESPONSE", "ESTROGEN_RESPONSE_EARLY", "HYPOXIA")),
    GeneRatio = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)))

################################################################################
# # plotting - Figure 3B
################################################################################
B <- ggplot(em_BV_B_short, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = log10padj)) +
  scale_color_gradientn(
    colours = c("blue", "white", "red", "red3"),
    values = scales::rescale(c(0, -log10(0.05), 3, 5)),
    limits = c(0, 5),
    name = "-log10(p_adj)"
  ) +
  scale_size_continuous(
    name = "Count",
    limits = c(8, 30),   # Set count range (min to max)
    range = c(2, 6)
  ) +
  scale_x_continuous(limits = c(0.05, 0.15), expand = c(0, 0)) +
  theme_bw() +
  labs(x = "GeneRatio", y = NULL, size = "Count") +
  theme(
    panel.grid = element_blank()
  )

P <- ggplot(em_PV_P_short, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = log10padj)) +
  scale_color_gradientn(
    colours = c("blue", "white", "red", "red3"),
    values = scales::rescale(c(0, -log10(0.05), 3, 5)),
    limits = c(0, 5),
    name = "-log10(p_adj)"
  ) +
  scale_size_continuous(
    name = "Count",
    limits = c(8, 30),   # Set count range (min to max)
    range = c(2, 6)
  ) +
  scale_x_continuous(limits = c(0.05, 0.15), expand = c(0, 0)) +
  theme_bw() +
  labs(x = "GeneRatio", y = NULL, size = "Count") +
  theme(
    panel.grid = element_blank()
  )

################################################################################
# # plotting - volcano plot in Figure 3C
################################################################################
log2FC_threshold <- 1
pvalue_threshold <- 0.05
BAP1vsPBRM1_all_markers$X <- rownames(BAP1vsPBRM1_all_markers)
BAP1vsPBRM1_all_markers$p_val_adj <- ifelse(BAP1vsPBRM1_all_markers$p_val_adj == 0, 1e-300, BAP1vsPBRM1_all_markers$p_val_adj)
BAP1vsPBRM1_all_markers$minusLog10Pvalue <- -log10(BAP1vsPBRM1_all_markers$p_val_adj)
gene_to_label <- BAP1vsPBRM1_all_markers[BAP1vsPBRM1_all_markers$X %in% c("VEGFA", "GRM8", "HAVCR1", "ANGPTL4", "GRIN2A", "HIF1A", "SOX6", "TGFBI", "SERPINE1", "PBX3", "FN1", "EPO", "GEM", "TMEM45A", "EPAS1", "EPAS1", "PGF", "LOX", "CDK6", "COL8A1", "PBX1", "CFH", "EPHA6"),]

volcano_plot <- ggplot(BAP1vsPBRM1_all_markers, aes(x = avg_log2FC, y = minusLog10Pvalue)) +
  geom_point(aes(color = ifelse(p_val_adj > pvalue_threshold | abs(avg_log2FC) < log2FC_threshold, 'Non-significant', ifelse(p_val_adj < pvalue_threshold & avg_log2FC > log2FC_threshold & pct.1 < 0.2, 'B_pct<0.2', ifelse(p_val_adj < pvalue_threshold & avg_log2FC < -log2FC_threshold & pct.2 < 0.2, 'P_pct<0.2', ifelse(p_val_adj < pvalue_threshold & avg_log2FC > log2FC_threshold & pct.1 > 0.2, "BAP1_up", 'PBRM1_up')))))) +
  scale_color_manual(values = c('BAP1_up' = "#F37E78", "PBRM1_up" = "#4189CB", 'B_pct<0.2' = "grey", 'P_pct<0.2' = "grey", 'Non-significant' = "grey90")) +
  geom_text_repel(
    data = gene_to_label, aes(label = X),
    fontface = "bold.italic",
    box.padding = 0.5,  # Increase spacing
    point.padding = 0.3,
    segment.color = 'black',  # Clear linking lines
    min.segment.length = 0,  # Allow very short segments
    segment.size = 0.4,  # Line thickness
    max.overlaps = 15,  # Allow more labels
    force = 1.2,  # Increase repulsion
    direction = "both"  # Allow flexible movement
  ) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed") +
  theme_minimal() +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Significant") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")
volcano_plot
ggsave("./output/figures/Volcano_plot_BAP1vsPBRM1_all_markers.svg", plot = volcano_plot, device = "svg", width = 6, height = 6, dpi = 300)

################################################################################
# # plotting - violin plot in Figure 3D and all signatures
################################################################################
df <- read.csv("./data/scRNA/metadata_PB.csv")
df <- read.csv("./data/scRNA/metadata_PB_hypoxia.csv")
df <- df[, 24:25]
df <- df[df$Genotype_collapse %in% c("BAP1_any", "PBRM1_any"),]

df$Genotype_collapse <- factor(df$Genotype_collapse, levels = c("PBRM1_any", "BAP1_any"))
# plot individual score of interest
p <- ggplot(df, aes(x=Genotype_collapse, y=signature_1HIF_meta_gene, fill=Genotype_collapse)) + 
  geom_violin() + 
  scale_fill_manual(values = c("#4189CB", "#F37E78")) +
  geom_boxplot(width=0.1) + stat_compare_means(method = "wilcox.test", comparisons = list(c("BAP1_any", "PBRM1_any")), label = "p.signif") + # Add pairwise comparisons p-value
  theme_minimal() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )
p
ggsave("./output/figures/Violin_HIF.svg", plot = p, device = "svg", width = 6, height = 8, dpi = 300)

features <- setdiff(names(df), "Genotype_collapse")[24:71]  # All columns except the grouping variable

# Create and save a violin plot for each feature
for (feature in features) {
  p <- ggplot(df, aes_string(x = "Genotype_collapse", y = feature, fill = "Genotype_collapse")) +
    geom_violin() +
    scale_fill_manual(values = c("#4189CB", "#F37E78")) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    stat_compare_means(
      method = "wilcox.test",
      comparisons = list(c("BAP1_any", "PBRM1_any")),
      label = "p.signif"
    ) +
    theme_minimal()
  
  # Export to PDF
  ggsave(filename = paste0("./violins/violin_", feature, ".pdf"), plot = p, width = 4, height = 4)
}


p <- ggplot(df, aes(x=Genotype_collapse, y=MP6_Hypoxia1, fill=Genotype_collapse)) + 
  geom_violin() + 
  scale_fill_manual(values = c("#4189CB", "#F37E78")) +
  geom_boxplot(width=0.1) + stat_compare_means(method = "wilcox.test", comparisons = list(c("BAP1_any", "PBRM1_any")), label = "p.signif")

################################################################################
# # plotting - fibroblast siganture heatmap in Figure 3H
################################################################################
meta <- read.csv(file = "./data/scRNA/fib_metadata.csv",row.names = 1)
geneset_totest <- read_csv("./data/scRNA/fib_marker.csv")
sig_list <- paste0("signature_1", colnames(geneset_totest))
# plot density and output mirror violin plots
for (sig in sig_list) {
  max_density <- max(density(meta[meta$Genotype_collapse == "BAP1_any", ][,sig])$y, 
                     density(meta[meta$Genotype_collapse == "PBRM1_any", ][, sig])$y) + 1
  # Plot
  p <- ggplot(meta, aes_string(x = sig, fill = "Genotype_collapse")) +
    geom_density(data = subset(meta, Genotype_collapse == "BAP1_any"),
                 aes(y = after_stat(density)), alpha = 0.6) +
    geom_density(data = subset(meta, Genotype_collapse == "PBRM1_any"),
                 aes(y = -after_stat(density)), alpha = 0.6) +
    scale_fill_manual(values = c("BAP1_any" = "#F37E78", "PBRM1_any" = "#4189CB")) +  # red and blue
    scale_x_continuous(limits = c(0, 0.9)) +  # 👈 set fixed x-axis
    scale_y_continuous(limits = c(-max_density, max_density), labels = abs) +
    labs(x = paste0("fib_", substr(sig, 12,14))) +
    geom_hline(yintercept = 0, color = "black") +         # Horizontal center line
    geom_vline(xintercept = 0, color = "black") +         # Vertical line at x = 0
    theme_minimal() +
    theme(
      panel.grid = element_blank(),   # removes both major and minor grid lines
      legend.position = "none"
    )
  ggsave(filename = paste0("./output/figures/fib_violins/violin_", sig, ".pdf"), plot = p, width = 6, height = 4)
}

# Calculate means for each group and generate df
meta1 <- meta[, 24:44]
grouped_means <- meta1 %>%
  group_by(Genotype_collapse) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  as.data.frame()  # Optional, for base R output
rownames(grouped_means) <- grouped_means$Genotype_collapse
grouped_means$Genotype_collapse <- NULL
colnames(grouped_means) <- c("MYH11+ SMC",	"ADAMDEC1+ fib",	"COL15A1+ progenitor fib",	"LRRC15+ myofib",	"PI16+ progenitor fib",	"ADH1B+ alveolar fib",	"IL6+ inflammatory fib",	"RGS5+ fib",	"CTNNB1+ myofib",	"PRG4+ synovial lining fib",	"HOPX+ myofib",	"MSLN+ mesothelial cells",	"HGF+ fib",	"HSPA6+ stress-response fib",	"SOX6+ epithelial crypt fib",	"SFRP2+ myofib",	"STMN1+ proliferative fib",	"HHIP+ SMC",	"MMP1+ fib",	"CD74+ antigen-presenting fib")
# filter out specilised fibroblast types (non-kidney tissue-specific)
grouped_means <- grouped_means[, c(-10, -18, -2, -15, -7, -6, -12)]
p <- pheatmap(grouped_means, 
              show_rownames = T, 
              show_colnames = T,
              cluster_rows = F,
              cluster_cols = T,
              color = colorRampPalette(c("#113c8a", "white", "#a84e0d"))(100),
              fontsize = 7
)
ggsave(filename = "./output/figures/fib_heatmap_map.pdf", p, width = 6, height = 3)
