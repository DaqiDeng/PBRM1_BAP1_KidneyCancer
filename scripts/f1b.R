library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
# working directory <- "PBRM1_BAP1_KidneyCancer"
setwd("path/to/PBRM1_BAP1_KidneyCancer")
################################################################################
# import data.frame and prepare clean data for plotting
oncodata <- read.csv("./data/fig1b_raw.csv", row.names=1)
df <- t(oncodata[, 7:32])
meta <- oncodata[, 1:6]
meta <- meta[, c(6:3)]

ann_colors = list(
  Patient = c(K1013 = "#E41A1C", K1022 = "#377EB8", K1030 = "#4DAF4A", K1040 = "#984EA3", K1042 = "#FF7F00", K1043 = "#DD9F28",  K1063 = "#FFFF33", K1184 = "#A65628"),
  Sample_Type = c(PDO = "wheat2", PDCL = "wheat4"),
  Tumour_Grade = c(G3 = "#f8992d", G4 = "#cc4e26"), # G1 = "#fcf9d6", G2 = "#ffd98e", 
  Tumour_Stage = c(T1a = "#bdd7e6", T3a = "#3082be", T4 = "#07519c"),
  Tumour_Necrosis = c(No = "lightgrey", Yes = "black"),
  Long_term_expansion = c(No = "lightgrey", Yes = "black")
)

df <- df[which(rowSums(df) > 0),]
df <- df[, c("K1013T2_GB1", "K1043T1_GB1", "K1022T1_GB1", "K1030T1_GB1", "K1030T2_GB1", "K1030T3_GB1", "K1063T1_GB1", "K1063T2_GB1", "K1040T1_GB1", "K1042T1_GB1", "K1184T1_GA1")]
# alteration color: 1 -> driver mutation, 2 -> CN loss, 3 -> CN gain 
A <- pheatmap(df,
         color = c("white", "#248444", "#0c92cf", "#ee4035"), # Colors for mutations
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = meta,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         gaps_col = c(1, 2, 3, 6, 8, 9, 10, 11),
         cellwidth = 12,  # Adjust cell width to increase spacing
         cellheight = 8,  # Adjust cell height to increase spacing
         fontsize = 6,
         fontsize_row = 6,
         fontsize_col = 6,
         border_color = "white",
         margins = c(8, 8))
# Save the heatmap to a file
ggsave("./output/figures/fig1b_heatmap.pdf", plot = A, width = 4, height = 3)
################################################################################
