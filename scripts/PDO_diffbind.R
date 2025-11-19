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
library(svglite)
################################################################################
# Define sample metadata and file locations
################################################################################
setwd("~/Downloads/PBRM1_BAP1_KidneyCancer/")
# the downloaded raw data including bam and peak files should be stored in the following structure, note this is using replicate-merged files
peak_dir <- "./data/PDO_TF/PDOB/peaks/"
peak_files <- list.files(peak_dir, pattern = "\\.narrowPeak", full.names = TRUE)

bam_dir <- "./data/PDO_TF/PDOB/bam/"
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)

sample_names <- c("K1013T2_GB1", "K1022N1_GN1", "K1022T1_GB1", "K1030N1_GN1", "K1030T1_GB1", "K1030T2_GB1", "K1030T3_GB1",  "K1040N1_GN1", "K1040T1_GB1", "K1063N1_GN1", "K1063T1_GB1", "K1063T2_GB1", "K1065N1_GN1", "K1184T1_GA1")
treatment <- c("V", "N", "VBP", "N", "VB", "VB", "VB", "N", "VP", "N", "VB", "VB", "N", "VP")

# Create the sample sheet data frame
sample_sheet <- data.frame(
  SampleID = sample_names,       # Unique identifier for each sample
  Condition = treatment,         # Experimental condition/group
  bamReads = bam_files,          # Path to aligned reads (BAM)
  Peaks = peak_files,            # Path to called peaks
  PeakCaller = "narrow"          # Type of peak files (narrow or broad)
)
# check the meta df is correct before constructing dna object, this can be saved as well.
################################################################################
# Create the DBA object from the sample sheet
################################################################################

atac_db <- dba(sampleSheet = sample_sheet)
dba.show(atac_db)

# Count reads overlapping peaks *!*computationally heavy step*!*
atac_db <- dba.count(atac_db,
                     summits = 250,       # Center peaks on summits, extend +/- 250bp
                     minOverlap = 2,      # Require at least 2 samples to share a peak
                     bUseSummarizeOverlaps = TRUE,
                     bParallel = TRUE)  # Use GenomicRanges::summarizeOverlaps()

dba.show(atac_db)

# Extract the raw count matrix (optional)
counts <- dba.peakset(atac_db, bRetrieve = TRUE)
head(counts)

# Generate correlation heatmap
dba.plotHeatmap(atac_db, correlations=F, maxSites = 5000,
                bUsePval=FALSE, th=1e-3, bRetrieve=TRUE)

# Create PCA plot to visualize sample relationships
dba.plotPCA(atac_db, DBA_CONDITION, label=DBA_ID, b3D = T)


mat_raw <- dba.plotHeatmap(atac_db, correlations=FALSE, maxSites = 5000,
                           bUsePval=FALSE, th=1e-3, bRetrieve=TRUE)
mat_raw <- as.data.frame(counts)
mat_raw <- as.data.frame(mat_raw)[, 6:19]

# compute sample correlation on raw (or variance-stabilised) values
#    Try Pearson (or spearman pearson if you expect outliers)
C <- cor(as.data.frame(mat_raw), method = "spearman", use = "pairwise.complete.obs")

# convert to a distance for clustering
D <- as.dist(1 - C)

# plot the correlation matrix but cluster using D
p <- pheatmap(
  C[c("K1040N1_GN1", "K1063N1_GN1", "K1022N1_GN1", "K1030N1_GN1", "K1065N1_GN1", "K1013T2_GB1", "K1184T1_GA1", "K1040T1_GB1", "K1022T1_GB1", "K1063T1_GB1", "K1063T2_GB1", "K1030T1_GB1", "K1030T2_GB1", "K1030T3_GB1"), c("K1040N1_GN1", "K1063N1_GN1", "K1022N1_GN1", "K1030N1_GN1", "K1065N1_GN1", "K1013T2_GB1", "K1184T1_GA1", "K1040T1_GB1", "K1022T1_GB1", "K1063T1_GB1", "K1063T2_GB1", "K1030T1_GB1", "K1030T2_GB1", "K1030T3_GB1")],
  color = colorRampPalette(c("white","#99d594","#1a9850"))(100),
  cluster_rows = F,
  cluster_cols = F,
  clustering_distance_rows = D,
  clustering_distance_cols = D,
  clustering_method = "ward.D2",   # try "average"; also test "ward.D2"
  display_numbers = FALSE,
  legend = TRUE,
  main = "Sample correlations"
)
ggsave("cor.svg", p, device = "svg",  width = 8, height = 6, dpi = 300)

mat <- dba.plotHeatmap(
  atac_db,
  correlations = FALSE,  # site-by-sample heatmap
  maxSites = 1000,
  scale        = "row",
  bUsePval     = FALSE,
  th           = 1e-3,
  bRetrieve    = TRUE     # <-- key: returns the matrix
)
mat1 <- as.data.frame(mat)[, 6:19]
# colour function: centred at 0 for row-scaled Z-scores
col_fun <- colorRamp2(c(-2, 0, 2), c("#07519c", "#F7F7F7", "#cc4e26"))

# Colour palette (centre = 0)
heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# Plot
pheatmap(scale(t(mat1)),
         color              = heat_colors,
         show_rownames       = TRUE,
         show_colnames       = FALSE,
         cluster_rows        = F,
         cluster_cols        = T,
         border_color        = NA,
         fontsize_col        = 9,
         treeheight_row      = 30,
         treeheight_col      = 30,
         main                = "Peak scores",
         legend              = TRUE,
)
