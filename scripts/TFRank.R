library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
setwd("~/Downloads/PBRM1_BAP1_KidneyCancer/")
# Load data
df <- read.csv("./data/PDO_TF/TFBS.csv", check.names = F)
# PRODUCE A RANK PLOT USING K1030_GN1
################################################################################
genes_to_annotate <- c("JUN::JUNB", "JUNB", "PAX8", "HNF1A", "HNF1B", "ELK1", "ETS1")

K1030N1_GN1 <- df[order(df$K1030N1_GN1_mean_score, decreasing = TRUE), ]
K1030N1_GN1$rank <- 1:nrow(K1030N1_GN1)

K1030N1_GN1$label <- ifelse(K1030N1_GN1$name %in% genes_to_annotate, K1030N1_GN1$name, NA)

p1 <- ggplot(K1030N1_GN1, aes(x = rank, y = K1030N1_GN1_mean_score)) +
  geom_line(color = "grey80") +
  geom_point(size = 0.7, alpha = 0.8) +
  geom_text_repel(
    aes(label = label),
    size = 6,
    max.overlaps = 10,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "red",       # line color
    segment.size = 1,              # line thickness
    min.segment.length = 0         # always draw line
  ) +
  labs(x = "TF Rank",
       y = "Binding Score") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  )
ggsave("./output/figures/TFBS_rank_K1030N1_GN1.svg", p1, device = "svg", width = 4, height = 6, dpi = 300)

K1013T2_GB1 <- df[order(df$K1013T2_GB1_mean_score, decreasing = TRUE), ]
K1013T2_GB1$rank <- 1:nrow(K1013T2_GB1)

K1013T2_GB1$label <- ifelse(K1013T2_GB1$name %in% genes_to_annotate, K1013T2_GB1$name, NA)

p2 <- ggplot(K1013T2_GB1, aes(x = rank, y = K1013T2_GB1_mean_score)) +
  geom_line(color = "grey80") +
  geom_point(size = 0.7, alpha = 0.8) +
  geom_text_repel(
    aes(label = label),
    size = 6,
    max.overlaps = 10,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "red",       # line color
    segment.size = 1,              # line thickness
    min.segment.length = 0         # always draw line
  ) +
  labs(x = "TF Rank",
       y = "Binding Score") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  )
ggsave("./output/figures/TFBS_rank_K1013T2_GB1.svg", p2, device = "svg", width = 4, height = 6, dpi = 300)

K1030T2_GB1 <- df[order(df$K1030T2_GB1_mean_score, decreasing = TRUE), ]
K1030T2_GB1$rank <- 1:nrow(K1030T2_GB1)

K1030T2_GB1$label <- ifelse(K1030T2_GB1$name %in% genes_to_annotate, K1030T2_GB1$name, NA)

p3 <- ggplot(K1030T2_GB1, aes(x = rank, y = K1030T2_GB1_mean_score)) +
  geom_line(color = "grey80") +
  geom_point(size = 0.7, alpha = 0.8) +
  geom_text_repel(
    aes(label = label),
    size = 6,
    max.overlaps = 10,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "red",       # line color
    segment.size = 1,              # line thickness
    min.segment.length = 0         # always draw line
  ) +
  labs(x = "TF Rank",
       y = "Binding Score") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  )
ggsave("./output/figures/TFBS_rank_K1030T2_GB1.svg", p3, device = "svg", width = 4, height = 6, dpi = 300)

K1040T1_GB1 <- df[order(df$K1040T1_GB1_mean_score, decreasing = TRUE), ]
K1040T1_GB1$rank <- 1:nrow(K1040T1_GB1)

K1040T1_GB1$label <- ifelse(K1040T1_GB1$name %in% genes_to_annotate, K1040T1_GB1$name, NA)

p4 <- ggplot(K1040T1_GB1, aes(x = rank, y = K1040T1_GB1_mean_score)) +
  geom_line(color = "grey80") +
  geom_point(size = 0.7, alpha = 0.8) +
  geom_text_repel(
    aes(label = label),
    size = 6,
    max.overlaps = 10,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "red",       # line color
    segment.size = 1,              # line thickness
    min.segment.length = 0         # always draw line
  ) +
  labs(x = "TF Rank",
       y = "Binding Score") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # remove major gridlines
    panel.grid.minor = element_blank(),   # remove minor gridlines
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # add black frame
    axis.line = element_line(color = "black") # add axes
  )
ggsave("./output/figures/TFBS_rank_K1040T1_GB1.svg", p4, device = "svg", width = 4, height = 6, dpi = 300)
