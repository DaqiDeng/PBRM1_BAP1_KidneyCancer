library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(GSVA)
library(GSVAdata)
library(SummarizedExperiment)
library(pheatmap)
library(RColorBrewer)
library(PCAtools)
library(ConsensusTME)
setwd("~/Downloads/PBRM1_BAP1_KidneyCancer/")
################################################################################
#TCGA data curation
################################################################################
# load TCGA_meta and clean   
TCGA_meta <- read_csv("./data/Bulk_RNA/TCGA_KIRC/all_meta_TCGA_final.csv")
TCGA_meta <- as.data.frame(TCGA_meta)
rownames(TCGA_meta) <- TCGA_meta$Sample
TCGA_meta <- TCGA_meta[, -1][,-1]
TCGA_meta <- TCGA_meta[TCGA_meta$SampleType == "Primary", ]
TCGA_meta <- TCGA_meta[order(row.names(TCGA_meta)), ]

# load TCGA_data count matrix and clean
TCGA_data <- read_csv("./data/Bulk_RNA/TCGA_KIRC/cts_FPKM_TCGA_final.csv")
TCGA_data <- as.data.frame(TCGA_data)
rownames(TCGA_data) <- TCGA_data$...1
TCGA_data <- TCGA_data[,-1]
TCGA_data <- TCGA_data[rownames(TCGA_meta), ]
TCGA_data <- TCGA_data[order(row.names(TCGA_data)), ]
################################################################################

################################################################################
#Tumourgraft data curation
################################################################################
# load TG_meta and clean   
TG_meta <- read_csv("./data/Bulk_RNA/JB_Tumourgraft/JB_TG_Cohort_VBP_meta.csv")
TG_meta <- as.data.frame(TG_meta)
rownames(TG_meta) <- TG_meta$RNAseq
TG_meta <- TG_meta[,0:21]
TG_meta <- TG_meta[order(row.names(TG_meta)), ]

# load TG_data count matrix and clean
TG_data <- read_csv("./data/Bulk_RNA/JB_Tumourgraft/JB_TG_Cohort_tpm.csv")
TG_data <- as.data.frame(TG_data)
TG_data <- TG_data[!is.na(TG_data$...1), ]
rownames(TG_data) <- TG_data$...1
TG_data <- TG_data[,-1]
TG_data <- as.data.frame(t(TG_data))
TG_data <- TG_data[rownames(TG_meta), ]
TG_data <- TG_data[order(row.names(TG_data)), ]
################################################################################

################################################################################
#TRACERx data curation
################################################################################
# load meta and clean
Tx_meta <- read_csv("./data/Bulk_RNA/TRACERxRenal/metadata_std.csv")
Tx_meta <- as.data.frame(Tx_meta)
rownames(Tx_meta) <- Tx_meta$Sample
Tx_meta <- Tx_meta[, -1]
Tx_meta <- Tx_meta[Tx_meta$SampleType != "Normal",]
Tx_meta <- Tx_meta[!is.na(Tx_meta$Grade),]
Tx_meta <- Tx_meta[order(row.names(Tx_meta)), ]

# load data count matrix and clean
Tx_data <- readRDS("./data/Bulk_RNA/TRACERxRenal/Tx_RNA_DEseq_vsd.rds")
Tx_data <- Tx_data@assays@data@listData[[1]]
Tx_data <- as.data.frame(t(Tx_data))
Tx_data <- Tx_data[rownames(Tx_meta),]
Tx_data <- Tx_data[order(row.names(Tx_data)), ]
################################################################################

################################################################################
#Filter data and subset for intersection
################################################################################
# subset intersect genes 
protein_coding_genes <- read.csv("./data/Bulk_RNA/protein_coding_genes.csv", sep="")
protein_coding_genes <- protein_coding_genes$symbol
intersect_genes <- intersect(colnames(TCGA_data), colnames(TG_data))
intersect_genes <- intersect(intersect_genes, colnames(Tx_data))
intersect_genes <- intersect(intersect_genes, protein_coding_genes)
TCGA_data <- TCGA_data[, intersect_genes]
TG_data <- TG_data[, intersect_genes]
Tx_data <- Tx_data[, intersect_genes]
rm(protein_coding_genes, intersect_genes)

# subset VHL, PBRM1, BAP1 COHORT
TCGA_meta <- TCGA_meta[TCGA_meta$VHL == 1 & (TCGA_meta$PBRM1 == 1 | TCGA_meta$BAP1 == 1), ]
TG_meta <- TG_meta[TG_meta$VHL == 1 & (TG_meta$PBRM1 == 1 | TG_meta$BAP1 == 1), ]
Tx_meta <- Tx_meta[Tx_meta$VHL == 1 & (Tx_meta$PBRM1 == 1 | Tx_meta$BAP1 == 1), ]
TCGA_data <- TCGA_data[rownames(TCGA_meta), ]
TG_data <- TG_data[rownames(TG_meta), ]
Tx_data <- Tx_data[rownames(Tx_meta), ]
all(rownames(TG_meta) == rownames(TG_data)) # should be TRUE
all(rownames(TCGA_meta) == rownames(TCGA_data)) # should be TRUE
all(rownames(Tx_meta) == rownames(Tx_data)) # should be TRUE
################################################################################
################################################################################
# Run ConsensusTME, TG cohort has high tumour purity (>90%) therefore no need to run
################################################################################
TCGA_TME <- ConsensusTME::consensusTMEAnalysis(t(TCGA_data), cancer = "KIRC", statMethod = "ssgsea")
Tx_TME <- ConsensusTME::consensusTMEAnalysis(t(Tx_data), cancer = "KIRC", statMethod = "ssgsea")

# Transpose so each row = sample, each column = cell type
TCGA_TME <- as.data.frame(t(TCGA_TME))
TCGA_TME$SampleID <- rownames(TCGA_TME)
Tx_TME <- as.data.frame(t(Tx_TME))
Tx_TME$SampleID <- rownames(Tx_TME)

write.csv(TCGA_TME, "./output/csv/TCGA_TME_scores.csv", row.names = FALSE)
write.csv(Tx_TME, "./output/csv/Tx_TME_scores.csv", row.names = FALSE)
# add VBP annotation
TCGA_meta_tme <- TCGA_meta %>%
  mutate(Group = case_when(
    VHL == 1 & BAP1 == 1 & PBRM1 == 0 ~ "VB",
    VHL == 1 & BAP1 == 0 & PBRM1 == 1 ~ "VP",
    VHL == 1 & BAP1 == 1 & PBRM1 == 1 ~ "VBP"
  ))
Tx_meta_tme <- Tx_meta %>%
  mutate(Group = case_when(
    VHL == 1 & BAP1 == 1 & PBRM1 == 0 ~ "VB",
    VHL == 1 & BAP1 == 0 & PBRM1 == 1 ~ "VP",
    VHL == 1 & BAP1 == 1 & PBRM1 == 1 ~ "VBP"
  ))

all(rownames(TCGA_TME) == rownames(TCGA_meta_tme))
all(rownames(Tx_TME) == rownames(Tx_meta_tme))

################################################################################
# heatmap in Figure S5C
################################################################################
# get group-wise scores
df_scores  <- t(TCGA_TME[,c(-19,-20)])   
sample_anno <- TCGA_meta 

P_df_score <- df_scores[, rownames(sample_anno[sample_anno$VHL == "1" & sample_anno$PBRM1 == "1" & sample_anno$BAP1 != "1",])]
B_df_score <- df_scores[, rownames(sample_anno[sample_anno$VHL == "1" & sample_anno$PBRM1 != "1" & sample_anno$BAP1 == "1",])]
PB_df_score <- df_scores[, rownames(sample_anno[sample_anno$VHL == "1" & sample_anno$PBRM1 == "1" & sample_anno$BAP1 == "1",])]

P <- pheatmap(P_df_score,
              cluster_rows = F, cluster_cols = T, 
              clustering_distance_rows = 'euclidean',
              clustering_method = 'ward.D',
              show_colnames = F)

B <- pheatmap(B_df_score,
              cluster_rows = F, cluster_cols = T, 
              clustering_distance_rows = 'euclidean',
              clustering_method = 'ward.D',
              show_colnames = F)

PB <- pheatmap(PB_df_score,
               cluster_rows = F, cluster_cols = T, 
               clustering_distance_rows = 'euclidean',
               clustering_method = 'ward.D',
               show_colnames = F)

P_df_score <- P_df_score[, P$tree_col[["order"]]]
B_df_score <- B_df_score[, B$tree_col[["order"]]]
PB_df_score <- PB_df_score[, PB$tree_col[["order"]]]

test_TCGA <- cbind(P_df_score, B_df_score, PB_df_score)
ann_colors = list(
  Purity = brewer.pal(9, "Blues"),
  VHL = c("1" = "#248444", "0" = "white"),
  PBRM1 = c("1" = "#248444", "0" = "white"),
  BAP1 = c("1" = "#248444", "0" = "white")
)
sample_anno$VHL <- as.factor(sample_anno$VHL)
sample_anno$PBRM1 <- as.factor(sample_anno$PBRM1)
sample_anno$BAP1 <- as.factor(sample_anno$BAP1)

pheatmap(test_TCGA,
         cluster_rows = T, cluster_cols = F, 
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = sample_anno[,c(9,11:13)],
         annotation_colors = ann_colors,
         show_colnames = F,
)

################################################################################
################################################################################
# Comparison TME cells composition - boxplots in Figure 3E S5D
################################################################################
# Add annotation (P/B)
TCGA_TME$PB <- TCGA_meta_tme$Group
Tx_TME$PB <- Tx_meta_tme$Group

TCGA_TME <- TCGA_TME[TCGA_TME$PB != "VBP", ]
TCGA_TME$PB <- factor(TCGA_TME$PB, levels = c("VP", "VB"))
Tx_TME <- Tx_TME[Tx_TME$PB != "VBP", ]
Tx_TME$PB <- factor(Tx_TME$PB, levels = c("VP", "VB"))

p <- ggplot(TCGA_TME, aes(x = PB, y = Fibroblasts, fill = PB)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("VP" = "#4189CB", "VB" = "#F37E78")) +
  labs(
    x = NULL,
    y = "Fibroblasts score"
  ) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("VB", "VP")), label = "p.format") + # Add pairwise comparisons p-values
  theme_classic(base_size = 14)
ggsave(filename = "./output/figures/TCGA_fib.svg", p, device = "svg", width = 4, height = 6, dpi = 300)


p <- ggplot(Tx_TME, aes(x = PB, y = Fibroblasts, fill = PB)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("VP" = "#4189CB", "VB" = "#F37E78")) +
  labs(
    x = NULL,
    y = "Fibroblasts score"
  ) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("VB", "VP")), label = "p.format") + # Add pairwise comparisons p-values
  theme_classic(base_size = 14)
ggsave(filename = "./output/figures/Tx_fib.svg", p, device = "svg", width = 4, height = 6, dpi = 300)

p <- ggplot(TCGA_TME, aes(x = PB, y = Endothelial, fill = PB)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("VP" = "#4189CB", "VB" = "#F37E78")) +
  labs(
    x = NULL,
    y = "Endothelial cell score"
  ) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("VB", "VP")), label = "p.format") + # Add pairwise comparisons p-values
  theme_classic(base_size = 14)
ggsave(filename = "./output/figures/TCGA_endo.pdf", p, device = "pdf", width = 5, height = 6, dpi = 300)

p <- ggplot(Tx_TME, aes(x = PB, y = Endothelial, fill = PB)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  scale_fill_manual(values = c("VP" = "#4189CB", "VB" = "#F37E78")) +
  labs(
    x = NULL,
    y = "Endothelial cell score"
  ) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("VB", "VP")), label = "p.signif") + # Add pairwise comparisons p-values
  theme_classic(base_size = 14)
ggsave(filename = "./output/figures/Tx_endo.pdf", p, device = "pdf", width = 5, height = 6, dpi = 300)
