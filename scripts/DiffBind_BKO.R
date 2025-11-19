library(data.table)
library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)
library(tidyverse)
library(ggrepel)
################################################################################
# Define sample metadata and file locations
################################################################################
setwd("~/Downloads/PBRM1_BAP1_KidneyCancer/")
# The downloaded raw data including bam and peak files should be in the right folder structure
peak_dir <- "./data/PDO_TF/BKO/peaks/"
peak_files <- list.files(peak_dir, full.names = TRUE)

bam_dir <- "./data/PDO_TF/BKO/bam/"
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

sample_names <- c("K1040T1_GB1_bko2_R1", "K1040T1_GB1_bko2_R2", "K1040T1_GB1_bko3_R1", "K1040T1_GB1_bko3_R2", "K1040T1_GB1_R1", "K1040T1_GB1_R2_1", "K1040T1_GB1_R2_2")
treatment <- c("bko", "bko", "bko", "bko", "bwt", "bwt", "bwt")

# Create the sample sheet data frame
sample_sheet <- data.frame(
  SampleID = sample_names,       # Unique identifier for each sample
  Condition = treatment,         # Experimental condition/group
  bamReads = bam_files,          # Path to aligned reads (BAM)
  Peaks = peak_files,            # Path to called peaks
  PeakCaller = "narrow"          # Type of peak files (narrow or broad)
)

# Save the sample sheet for future reference, check the meta is correct 
write.csv(sample_sheet, "./scripts/atac_seq_sample_sheet_bko.csv", row.names = FALSE)

################################################################################
# Create the DBA object from the sample sheet
################################################################################
atac_db <- dba(sampleSheet = sample_sheet)

# Examine the DBA object
dba.show(atac_db)

# Count reads overlapping peaks
atac_db <- dba.count(atac_db,
                     summits = 250,       # Center peaks on summits, extend +/- 250bp
                     minOverlap = 2,      # Require at least 2 samples to share a peak
                     bUseSummarizeOverlaps = TRUE)  # Use GenomicRanges::summarizeOverlaps()

# View updated DBA object with count information
dba.show(atac_db)

# Extract the raw count matrix (optional)
counts <- dba.peakset(atac_db, bRetrieve = TRUE)
head(counts)
# save the object 
saveRDS(atac_db, file = "./atac_db_bko.rds")
atac_db <- readRDS("./atac_db_bko.rds")

# Generate correlation heatmap
dba.plotHeatmap(atac_db)

################################################################################
# Differential analysis
################################################################################

# Define the contrast by condition
atac_db <- dba.contrast(atac_db, 
                        categories = DBA_CONDITION,  # Use the Condition column from sample sheet
                        minMembers = 2)              # Minimum samples per group

# Display the contrast information
dba.show(atac_db, bContrasts = TRUE)

# Run differential analysis (default: DESeq2)
atac_db <- dba.analyze(atac_db)

# Display analysis results summary
dba.show(atac_db, bContrasts = TRUE)

# Default: get significant results (FDR < 0.001)
diff_peaks <- dba.report(atac_db, th = 0.001)


# Convert to data frame and write to CSV
diff_df <- as.data.frame(diff_peaks)
fwrite(diff_df, "./output/csv/differential_peaks_FDR001_bko.csv")


# upregulated peaks
diff_peaks_up <- diff_peaks[diff_peaks$Fold > 0]
diff_df_up <- as.data.frame(diff_peaks_up)
# downregulated peaks
diff_peaks_down <- diff_peaks[diff_peaks$Fold < 0]
diff_df_down <- as.data.frame(diff_peaks_down)
################################################################################
# Annotate the differential peaks with nearby genes
################################################################################
# Use TxDb.Hsapiens.UCSC.hg38.knownGene and org.Hs.eg.db for human samples
diff_peaks_anno <- annotatePeak(diff_peaks,
                                TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,  # Mouse genome
                                annoDb = "org.Hs.eg.db",                     # Mouse gene annotations
                                tssRegion = c(-3000, 3000),                  # Define promoter region
                                verbose = FALSE)

# View detailed annotation table
annotated_df <- as.data.frame(diff_peaks_anno)
head(annotated_df)

# Save annotated results
fwrite(annotated_df, "./output/csv/differential_peaks_annotated_bko.csv")

################################################################################
# Visualizations Heatmap of differential peaks across samples
################################################################################
dba.plotHeatmap(atac_db, contrast=1, correlations=FALSE, scale="row", bUsePval=FALSE, th = 0.001)

mat <- dba.plotHeatmap(
  atac_db,
  contrast     = 1,
  correlations = FALSE,  # site-by-sample heatmap
  scale        = "row",
  bUsePval     = FALSE,
  th           = 1e-3,
  bRetrieve    = TRUE     # <-- key: returns the matrix
)
mat1 <- as.data.frame(mat)[, 6:12]

# Colour palette (centre = 0)
heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# Plot
pheatmap(t(scale(t(mat1))),
         color              = heat_colors,
         show_rownames       = FALSE,
         show_colnames       = TRUE,
         cluster_rows        = FALSE,
         cluster_cols        = FALSE,
         border_color        = NA,
         fontsize_col        = 9,
         treeheight_row      = 30,
         treeheight_col      = 30,
         main                = "Top1000 DARs",
         legend              = TRUE,
         )
################################################################################
# Visualizations - pie chart 
################################################################################
diff_peaks_anno_up <- annotatePeak(diff_peaks_up,
                                TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,  # Mouse genome
                                annoDb = "org.Hs.eg.db",                     # Mouse gene annotations
                                tssRegion = c(-3000, 3000),                  # Define promoter region
                                verbose = FALSE)
diff_peaks_anno_down <- annotatePeak(diff_peaks_down,
                                TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,  # Mouse genome
                                annoDb = "org.Hs.eg.db",                     # Mouse gene annotations
                                tssRegion = c(-3000, 3000),                  # Define promoter region
                                verbose = FALSE)

plotAnnoPie(diff_peaks_anno_up)
plotAnnoPie(diff_peaks_anno_down)

# collapse ChIPseeker annotation into broad classes
collapse_anno <- function(peakAnno) {
  df <- as.data.frame(peakAnno)
  df$Class <- case_when(
    grepl("^Promoter", df$annotation)                  ~ "Promoter",
    grepl("5' UTR", df$annotation)                    ~ "5' UTR",
    grepl("3' UTR", df$annotation)                    ~ "3' UTR",
    grepl("Exon",   df$annotation)                    ~ "Exon",
    grepl("Intron", df$annotation)                    ~ "Intron",
    grepl("^Downstream", df$annotation, ignore.case=TRUE) ~ "Downstream",
    grepl("Intergenic", df$annotation, ignore.case=TRUE)  ~ "Intergenic",
    TRUE ~ "Other"
  )
  # summarise counts + percentages
  df %>%
    dplyr::count(Class, name = "n") %>%
    mutate(pct = 100*n/sum(n)) %>%
    arrange(desc(n))
}

# helper: make a pie with ggplot (legend shown)
make_pie <- function(summary_df, title = NULL, palette = NULL) {
  # order factor so legend & pie are consistent
  summary_df$Class <- factor(summary_df$Class,
                             levels = c("Promoter","Exon","Intron","5' UTR","3' UTR",
                                        "Downstream","Intergenic","Other"))
  summary_df <- summary_df %>% filter(!is.na(Class))
  ggplot(summary_df, aes(x = "", y = n, fill = Class)) +
    geom_col(color = "black", width = 1) +
    coord_polar(theta = "y") +
    labs(title = title, x = NULL, y = NULL, fill = "Annotation") +
    scale_fill_manual(values = palette, drop = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

diff_peaks_anno_up1 <- as.data.frame(diff_peaks_anno_up)
diff_peaks_anno_down1 <- as.data.frame(diff_peaks_anno_down)
sum1 <- collapse_anno(diff_peaks_anno_up1)
sum2 <- collapse_anno(diff_peaks_anno_down1)

# consistent colors across charts (named vector; adjust if you like)
pal <- c(
  "Promoter"   = "#51ADF7",
  "Exon"       = "#E07B91",
  "Intron"     = "#C5ADD7",
  "5' UTR"     = "#7FC97F",
  "3' UTR"     = "#FBB4AE",
  "Downstream" = "#9E0142",
  "Intergenic" = "#FFF6A6",
  "Other"      = "#BDBDBD"
)

pal <- brewer.pal(8, "Set2")
names(pal) <- c("Promoter","Exon","Intron","5' UTR","3' UTR",
                "Downstream","Intergenic","Other")

p1 <- make_pie(sum1, title = "Increased accessible regions", palette = pal)
p2 <- make_pie(sum2, title = "Decreased accessible regions", palette = pal)

# stack vertically (up–down) with a single shared legend
legend <- cowplot::get_legend(p1 + theme(legend.position = "right"))
p1_noleg <- p1 + theme(legend.position = "none")
p2_noleg <- p2 + theme(legend.position = "none")

stacked <- cowplot::plot_grid(
  p1_noleg, p2_noleg, ncol = 1, align = "v", rel_heights = c(1, 1)
)

final_plot <- cowplot::plot_grid(stacked, legend, ncol = 2, rel_widths = c(1, 0.35))
final_plot

################################################################################
# Functional annotation of the peaks 
################################################################################
# Extract genes associated with differential peaks 
# (focusing on those within 5kb of a TSS)
diff_genes_up_promter <- unique(annotated_df$geneId[annotated_df$distanceToTSS <= 10000 & annotated_df$Fold > 0])
promoter_up <- annotated_df %>%
  filter(grepl("^Promoter", annotation),
         Fold > 0)
genes_promoter_up <- unique(promoter_up$SYMBOL)  # or $geneId if you prefer ENTREZ
gene_info <- AnnotationDbi::select(org.Hs.eg.db,
                                   keys = genes_promoter_up,
                                   keytype = "SYMBOL",
                                   columns = c("SYMBOL", "GENETYPE"))
protein_coding <- gene_info %>%
  filter(GENETYPE == "protein-coding") %>%
  pull(SYMBOL) %>%
  unique()

writeLines(as.character(na.omit(protein_coding)), "./data/PDO_TF/BKO/increased_accessibile_promoter_genes.txt")

enrichment_hallmark <- read.csv("./data/PDO_TF/BKO/enrichment_hallmark.csv")

enrichment_hallmark$Pathway <- substr(enrichment_hallmark$Pathway, 10, 50)
colnames(enrichment_hallmark) <- c("FDR", "nGenes", "Pathway_gene", "Fold_Enrichment", "Hallmark", "URL", "Genes")
df <- enrichment_hallmark %>%
  arrange(desc(Fold_Enrichment)) %>%
  mutate(Rank = row_number())
df <- df %>%
  mutate(Hallmark = fct_reorder(Hallmark, Fold_Enrichment, .desc = FALSE))  # highest first


p <- ggplot(df, aes(x = Fold_Enrichment,
                       y = Hallmark,
                       size = nGenes,
                       colour = FDR)) +
  geom_point(alpha = 0.95) +
  scale_colour_viridis_c(
    option = "D", 
    direction = -1,
    name = "FDR",
    trans = "log10",
    limits = c(min(df$FDR), max(df$FDR))
  ) +
  scale_size_area(max_size = 9, name = "Gene count") +
  labs(
    x = "Fold enrichment",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    plot.title   = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave("hallmark_enrichment.svg", p, device = "svg",  width = 8, height = 8, dpi = 300)

# --------------------------------------------
# Generate a volcano plot with the TFBS sites 
# --------------------------------------------
TFBS <- read.csv("./data/PDO_TF/BKO/bindetect_results_volcano.csv", check.names = FALSE)
# K1040_GE1_CTRL is K1040T1_GB1 (wildtype) and K1040_GE1 is K1040T1_GB1_bko (BAP1 knockout)
colnames(TFBS) <- c("TF", "K1040_GE1_CTRL_mean_score", "K1040_GE1_mean_score", "K1040_GE1_CTRL_K1040_GE1_change", "pvalue", "highlighted")
TFBS$DiffBS <- -TFBS$K1040_GE1_CTRL_K1040_GE1_change
TFBS$minusLog10Pvalue <- -log10(TFBS$pvalue)

label_TFs <- c(
  "FOSB::JUNB", "JUNB", "FOSL1::JUNB", "FOS::JUN", "BATF3", "FOS::JUND",
  "ZEB1", "NFIC", "NFIX", "CEBPB", "POU1F1", "POU3F2", "POU3F4", "POU5F1B", "POU2F3"
)

# Make sure types are correct
TFBS <- TFBS %>%
  mutate(
    highlighted = as.logical(highlighted),
    sign = if_else(DiffBS >= 0, "pos", "neg"),
    label_flag = TF %in% label_TFs
  )

p <- ggplot(TFBS, aes(x = DiffBS, y = minusLog10Pvalue)) +
  geom_point(
    aes(colour = sign, size = highlighted),
    alpha = 0.8, stroke = 0, shape = 16
  ) +
  geom_text_repel(
    data = filter(TFBS, label_flag),
    aes(label = TF, colour = sign),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.color = "grey30",      # connector line colour
    segment.size = 0.3,            # connector thickness
    segment.linetype = "solid",   # straight solid line
    min.segment.length = 0,       # ensure connectors always drawn
    force = 5,                     # stronger repulsion
    max.overlaps = Inf
  ) +
  scale_colour_manual(
    values = c(pos = "#F37E78", neg = "#4189CB"),
    breaks = c("pos","neg"),
    labels = c("High in BKO","High in WT"),
    name = "Direction"
  ) +
  scale_size_manual(
    values = c(`FALSE` = 1.8, `TRUE` = 3.8),
    breaks = c(FALSE, TRUE),
    labels = c("No","Yes"),
    name = "Significance"
  ) +
  labs(
    x = "Differential binding score",
    y = expression(-log[10](pvalue))
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

p
ggsave("TFBS_volcano_plot.svg", p, device = "svg",  width = 8, height = 8, dpi = 300)




