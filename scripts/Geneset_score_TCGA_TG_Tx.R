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
#Calculate ssGSEA scores
################################################################################
#construct a SummerizedExperiment object
TCGA_CTS <- SummarizedExperiment(assays = list(TPM = t(TCGA_data)), colData = TCGA_meta)
TG_CTS <- SummarizedExperiment(assays = list(TPM = t(TG_data)), colData = TG_meta)
Tx_CTS <- SummarizedExperiment(assays = list(TPM = t(Tx_data)), colData = Tx_meta)
# import the genesets 
geneset <- read_csv("./data/Bulk_RNA/genesets.csv")
new_geneset <- list()
for (i in 1:length(geneset)){
  new_geneset[[i]] <- geneset[[i]][!is.na(geneset[[i]])]
}
names(new_geneset) <- names(geneset)

rm(i, geneset)

# calculate GSEA Scores for the two datasets
TCGA_GSVA <- gsva(TCGA_CTS, new_geneset, method = "ssgsea")
TG_GSVA <- gsva(TG_CTS, new_geneset, method = "ssgsea")
Tx_GSVA <- gsva(Tx_CTS, new_geneset, method = "ssgsea")

# export sample annotations and clean up
TCGA_sample_anno_1 <- as.data.frame(TCGA_GSVA@colData) # get sample annotation
TCGA_sample_anno <- as.data.frame(lapply(TCGA_sample_anno_1, as.character))
rownames(TCGA_sample_anno) <- rownames(TCGA_sample_anno_1)
TCGA_sample_anno <- TCGA_sample_anno[, c(4:5, 9, 11:36)]
TCGA_sample_anno$Purity  <- as.numeric(TCGA_sample_anno$Purity)
rm(TCGA_sample_anno_1)

TG_sample_anno_1 <- as.data.frame(TG_GSVA@colData) # get sample annotation
TG_sample_anno <- as.data.frame(lapply(TG_sample_anno_1, as.character))
rownames(TG_sample_anno) <- rownames(TG_sample_anno_1)
TG_sample_anno <- TG_sample_anno[, c(6:10, 12, 17:21)]
TG_sample_anno$Purity  <- as.numeric(TG_sample_anno$Purity)
rm(TG_sample_anno_1)

Tx_sample_anno_1 <- as.data.frame(Tx_GSVA@colData) # get sample annotation
Tx_sample_anno <- as.data.frame(lapply(Tx_sample_anno_1, as.character))
rownames(Tx_sample_anno) <- rownames(Tx_sample_anno_1)
Tx_sample_anno <- Tx_sample_anno[, c(4:5, 9, 11:36)]
Tx_sample_anno$Purity  <- as.numeric(Tx_sample_anno$Purity)
rm(Tx_sample_anno_1)

# export ssGSEA scores
TCGA_scores <- as.data.frame(TCGA_GSVA@assays@data$es)
TG_scores <- as.data.frame(TG_GSVA@assays@data$es)
Tx_scores <- as.data.frame(Tx_GSVA@assays@data$es)

# scale each gsea score across the cohort
TCGA_scores <- as.data.frame(t(scale(t(TCGA_scores))))
TG_scores <- as.data.frame(t(scale(t(TG_scores))))
Tx_scores <- as.data.frame(t(scale(t(Tx_scores))))

################################################################################
#Export results
################################################################################
write.csv(TCGA_scores,"./output/csv/TCGA_ssGSEA.csv") # TCGA sample naming connected "-", note this in following analysis
write.csv(TCGA_sample_anno,"./output/csv/TCGA_ssGSEA_anno.csv") # TCGA sample naming connected "-", note this in following analysis
write.csv(TG_scores,"./output/csv/TG_ssGSEA.csv")
write.csv(TG_sample_anno,"./output/csv/TG_ssGSEA_anno.csv")
write.csv(Tx_scores,"./output/csv/Tx_ssGSEA.csv")
write.csv(Tx_sample_anno,"./output/csv/Tx_ssGSEA_anno.csv")

################################################################################


################################################################################
#Plot 
################################################################################
# get group-wise scores
df_scores  <- TG_scores   # change for three cohorts
sample_anno <- TG_sample_anno # change for three cohorts

P_df_score <- df_scores[, rownames(sample_anno[sample_anno$VHL == "1" & sample_anno$PBRM1 == "1" & sample_anno$BAP1 != "1",])]
B_df_score <- df_scores[, rownames(sample_anno[sample_anno$VHL == "1" & sample_anno$PBRM1 != "1" & sample_anno$BAP1 == "1",])]
PB_df_score <- df_scores[, rownames(sample_anno[sample_anno$VHL == "1" & sample_anno$PBRM1 == "1" & sample_anno$BAP1 == "1",])]

P <- pheatmap(P_df_score,
              cluster_rows = F, cluster_cols = T, 
              clustering_distance_rows = 'euclidean',
              clustering_method = 'ward.D',
              annotation_col = sample_anno[, c(1,2)],
              show_colnames = F)

B <- pheatmap(B_df_score,
         cluster_rows = F, cluster_cols = T, 
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = sample_anno[, c(1,2)],
         show_colnames = F)

PB <- pheatmap(PB_df_score,
              cluster_rows = F, cluster_cols = T, 
              clustering_distance_rows = 'euclidean',
              clustering_method = 'ward.D',
              annotation_col = sample_anno[, c(1,2)],
              show_colnames = F)

P_df_score <- P_df_score[, P$tree_col[["order"]]]
B_df_score <- B_df_score[, B$tree_col[["order"]]]
PB_df_score <- PB_df_score[, PB$tree_col[["order"]]]

test_TCGA <- cbind(P_df_score, B_df_score, PB_df_score) # change for three cohorts
test_TG <- cbind(P_df_score, B_df_score, PB_df_score) # change for three cohorts
test_Tx <- cbind(P_df_score, B_df_score, PB_df_score) # change for three cohorts

test <- test_TG
#combined metadata/sample annotation 
# sample_anno <- rbind(TCGA_sample_anno, Tx_sample_anno)
# TG_sample_anno[, setdiff(colnames(TCGA_sample_anno), colnames(TG_sample_anno))] <- NA
# sample_anno[, setdiff(colnames(TG_sample_anno), colnames(TCGA_sample_anno))] <- NA
# sample_anno$Source <- "Primary Tumor"
# sample_anno <- rbind(sample_anno, TG_sample_anno)
# sample_anno <- sample_anno[colnames(test),]

ann_colors = list(
  Grade = c("NA" = "lightgrey", GX = "lightgrey", G1 = "#fcf9d6", G2 = "#ffd98e", G3 = "#f8992d", G4 = "#cc4e26"),
  Stage = c("NA" = "lightgrey", T1 = "#bdd7e6", T1a = "#bdd7e6", T1b = "#bdd7e6", T1c = "#bdd7e6", T2 = "#6cadd7",T2a = "#6cadd7", T2b = "#6cadd7", T2c = "#6cadd7", T3 = "#3082be", T3a = "#3082be", T3b = "#3082be", T3c = "#3082be", T4 = "#07519c"),
  Sarcomatoid = c("NA" = "lightgrey", Present = "#dd2f28", Absent = "#fde0d2"),
  Rhabdoid = c("NA" = "lightgrey", Present = "#dd2f28", Absent = "#fde0d2"),
  # Source = c(Primary Tumor = "wheat2", Metastasis = "wheat4"),
  Purity = brewer.pal(9, "Blues"),
  VHL = c("1" = "#248444", "0" = "white"),
  PBRM1 = c("1" = "#248444", "0" = "white"),
  BAP1 = c("1" = "#248444", "0" = "white"),
  SETD2 = c("1" = "#248444", "0" = "white"),
  PIK3CA = c("1" = "#248444", "0" = "white"),
  MTOR = c("1" = "#248444", "0" = "white"),
  PTEN = c("1" = "#248444", "0" = "white"),
  KDM5C = c("1" = "#248444", "0" = "white"),
  TP53 = c("1" = "#248444", "0" = "white"),
  TSC1 = c("1" = "#248444", "0" = "white"),
  TSC2 = c("1" = "#248444", "0" = "white"),
  ARID1A = c("1" = "#248444", "0" = "white"),
  gain_1q = c("1" = "#ee4035", "0" = "white"),
  gain_2q = c("1" = "#ee4035", "0" = "white"),
  gain_5q = c("1" = "#ee4035", "0" = "white"),
  gain_7q = c("1" = "#ee4035", "0" = "white"),
  gain_8q = c("1" = "#ee4035", "0" = "white"),
  gain_12p = c("1" = "#ee4035", "0" = "white"),
  gain_20q = c("1" = "#ee4035", "0" = "white"),
  loss_1p = c("1" = "#0c92cf", "0" = "white"),
  loss_3p = c("1" = "#0c92cf", "0" = "white"),
  loss_4q = c("1" = "#0c92cf", "0" = "white"),
  loss_6q = c("1" = "#0c92cf", "0" = "white"),
  loss_8p = c("1" = "#0c92cf", "0" = "white"),
  loss_9p = c("1" = "#0c92cf", "0" = "white"),
  loss_14q = c("1" = "#0c92cf", "0" = "white")
)
################################################################################
# specifically for Tx cohort - heatmap Figure S3
################################################################################
sample_anno1 <- sample_anno[, c(-14, -30, -31)]
pheatmap(test,
         cluster_rows = F, cluster_cols = F, 
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = sample_anno1,
         annotation_colors = ann_colors,
         breaks = seq(-5, 5, length.out = 101),
         show_colnames = F,
)
################################################################################
# specifically for TG cohort - heatmap Figure S4
################################################################################
sample_anno1 <- sample_anno[, c(1:8, 11:10)]
pheatmap(test,
         cluster_rows = F, cluster_cols = F, 
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = sample_anno1,
         annotation_colors = ann_colors,
         breaks = seq(-3, 3, length.out = 101),
         show_colnames = F,
)

################################################################################
#Plot individual SCORES - FOR TCGA and Tx, with VBP
################################################################################
# plot individual score of interest
P_df_score["PB_status",] <- "VP"
B_df_score["PB_status",] <- "VB"
PB_df_score["PB_status",] <- "VPB"
df <- cbind(P_df_score, B_df_score, PB_df_score)
df <- t(df)
df <- as.data.frame(df)
df[, 1:42] <- as.numeric(unlist(df[,1:42]))

p <- ggplot(df, aes(x=PB_status, y=HALLMARK_TGF_BETA_SIGNALING, fill=PB_status)) + 
  geom_violin(alpha=0.5) +  
  geom_boxplot(width=0.2, outlier.shape = NA) + 
  scale_fill_manual(values=c("#F37E78", "#4189CB", "#8D1755")) + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("VB", "VP"), c("VP", "VPB"), c("VB", "VPB")),
                     label = "p.format",
                     label.y.npc = "top", 
                     aes(label = paste0("p = ", ..p.format..))) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  ) 
p
# loop for all signatures
# Define the output directory
output_dir <- "./output/figures/Boxplots/"

# Get the numeric columns (excluding the grouping variable 'PB_status')
numeric_cols <- setdiff(names(df), "PB_status")

# Loop through each column
for (col in numeric_cols) {
  
  p <- ggplot(df, aes_string(x = "PB_status", y = col, fill = "PB_status")) + 
    geom_violin(alpha = 0.5) +  
    geom_boxplot(width = 0.2, outlier.shape = NA) + 
    scale_fill_manual(values = c("#F37E78", "#4189CB", "#8D1755")) + 
    stat_compare_means(method = "wilcox.test", 
                       comparisons = list(c("VB", "VP"), c("VP", "VPB"), c("VB", "VPB")),
                       label = "p.format", 
                       label.y.npc = "top") +  
    theme_minimal()+
    theme(
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    ) 
  # Save the plot with a dynamic filename
  ggsave(filename = paste0(output_dir, paste0("Tx_", col), ".pdf"), 
         plot = p, device = "pdf", width = 4, height = 5, units = "in", dpi = 300)
}


################################################################################
#Plot individual SCORES - FOR TG, without VBP
################################################################################

### discard PB in TG cohort as only have 2 samples
P_df_score["PB_status",] <- "VP"
B_df_score["PB_status",] <- "VB"
PB_df_score["PB_status",] <- "VPB"
df <- cbind(P_df_score, B_df_score)
df <- t(df)
df <- as.data.frame(df)
df[, 1:43] <- as.numeric(unlist(df[,1:43]))
p <- ggplot(df, aes(x=PB_status, y=HALLMARK_TGF_BETA_SIGNALING, fill=PB_status)) +
  geom_violin() + 
  geom_boxplot(width=0.2, outlier.shape = NA) +
  scale_fill_manual(values=c("#F37E78", "#4189CB")) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("VB", "VP")), label = "p.format") + # Add pairwise comparisons p-values
  theme_minimal()+
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  ) 
p
# loop for all signatures
# Define the output directory
output_dir <- "./output/figures/Boxplots/"

# Get the numeric columns (excluding the grouping variable 'PB_status')
numeric_cols <- setdiff(names(df), "PB_status")

# Loop through each column
for (col in numeric_cols) {
  
  p <- ggplot(df, aes_string(x = "PB_status", y = col, fill = "PB_status")) + 
    geom_violin(alpha = 0.5) +  
    geom_boxplot(width = 0.2, outlier.shape = NA) + 
    scale_fill_manual(values = c("#F37E78", "#4189CB")) + 
    stat_compare_means(method = "wilcox.test", 
                       comparisons = list(c("VB", "VP")),
                       label = "p.format", 
                       label.y.npc = "top") +  
    theme_minimal()+
    theme(
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    ) 
  # Save the plot with a dynamic filename
  ggsave(filename = paste0(output_dir, paste0("TG_", col), ".pdf"), 
         plot = p, device = "pdf", width = 4, height = 5, units = "in", dpi = 300)
}

################################################################################
# for K143 only
################################################################################
K143_score <- Tx_scores[, grepl("^K153", names(Tx_scores))]
K143_sample_anno <- Tx_sample_anno[grepl("^K153", rownames(Tx_sample_anno)), ]

pheatmap(K143_score,
         cluster_rows = F, cluster_cols = T, 
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         annotation_col = K143_sample_anno,
         annotation_colors = ann_colors,
         show_colnames = T,
)

df <- K143_score
df["PB_status",] <- c("VB", "VB", "VP", "VP", "VP", "VB")
df <- t(df)
df <- as.data.frame(df)

# use tx 
DF_153 <- df[c("K153_R54", "K153_R55", "K153_R58", "K153_R26", "K153_R4", "K153_R59"),]

output_dir <- "./output/figures/Boxplots/"

# Get the numeric columns (excluding the grouping variable 'PB_status')
numeric_cols <- setdiff(names(DF_153), "PB_status")

# Loop through each column
for (col in numeric_cols) {
  
  p <- ggplot(DF_153, aes_string(x = "PB_status", y = col, fill = "PB_status")) + 
    geom_boxplot(width = 0.5) +
    scale_fill_manual(values = c("#F37E78", "#4189CB")) + 
    stat_compare_means(method = "wilcox.test", 
                       comparisons = list(c("VB", "VP")),
                       label = "p.format", 
                       label.y.npc = "top") +
    geom_dotplot(binaxis='y', stackdir='center') +
    theme_minimal()+
    theme(
      axis.text.x = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold")
    ) 
  
  # Save the plot with a dynamic filename
  ggsave(filename = paste0(output_dir, paste0("K153_", col), ".pdf"), 
         plot = p, device = "pdf", width = 4, height = 5, units = "in", dpi = 300)
}

