library(dplyr)
library(devtools)
library(curl)
library(tidyr)
library(Seurat)
library(SeuratObject)
library(scCustomize)
library(Nebulosa)
library(ggplot2)
library(purrr)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(UCell)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler) 
library(msigdbr)
library(fgsea)
library(enrichR) #need internet connection
library(Cairo)
# activate X11 display in Rstudio Server if on HPC
options(bitmapType='cairo')
setwd("~/Downloads/PBRM1_BAP1_KidneyCancer/")
########################################################################## 
# Read in objects from data folder
##########################################################################
PDOT <- readRDS("/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/Seurats/PDOT_v_deposit.rds")

PDOT <- SetIdent(PDOT, value = PDOT@meta.data$LineID)

PDOT_D4 <- DotPlot(PDOT, features = HIF_meta_gene)
ggsave(filename = "./PDOT_dotplot_hif.pdf", PDOT_D4, device = "pdf", width = 40, height = 8, dpi = 300)

PDOT_D5 <- DotPlot(PDOT, features = gene_of_interest)
ggsave(filename = "./PDOT_dotplot_GOI.pdf", PDOT_D5, device = "pdf", width = 32, height = 8, dpi = 300)


#markers by clusters - heatmap won't work
PDOT_Cluster_markers <- FindAllMarkers(PDOT, only.pos = TRUE, min.pct = 0.8, logfc.threshold = 0.25)
PDOT_Cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> PDOT_top10
PDOT_D3 <- DoHeatmap(PDOT, features = PDOT_top10$gene)

########################################################################## 
# Differential Expression Analysis
##########################################################################
## markers by EVO_TYPE
# N V VB VP VBP Specific
Idents(PDOT) <- "Evo_Type"
Makers <- FindAllMarkers(PDOT, min.pct = 0.25)
write.csv(Makers, file = "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/DEGs/Find_all_markers_Evo_Type.csv")

# BAP1 and PBRM1 compared to normal
Idents(PDOT) <- "Evo_Type"
BAP1vN_markers <- FindMarkers(PDOT, ident.1 = "BAP1_Driven", ident.2 = "Normal", max.cells.per.ident = 500, min.pct = 0.25, logfc.threshold = 0.5)
BAP1vN_markers <- BAP1vN_markers[BAP1vN_markers$p_val_adj < 0.05,]
write.csv(BAP1vN_markers, file = "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/DEGs/BAP1vsN_all_markers.csv")

PBRM1vN_markers <- FindMarkers(PDOT, ident.1 = "PBRM1_Driven", ident.2 = "Normal", max.cells.per.ident = 500, min.pct = 0.25, logfc.threshold = 0.5)
PBRM1vN_markers <- PBRM1vN_markers[PBRM1vN_markers$p_val_adj < 0.05,]
write.csv(PBRM1vN_markers, file = "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/DEGs/PBRM1vsN_all_markers.csv")

# BAP1 and PBRM1 compared to VHL monodriver
BAP1vV_markers <- FindMarkers(PDOT, ident.1 = "BAP1_Driven", ident.2 = "VHL_monodriver", max.cells.per.ident = 500, min.pct = 0.25, logfc.threshold = 0.5)
BAP1vV_markers <- BAP1vV_markers[BAP1vV_markers$p_val_adj < 0.05,]
write.csv(BAP1vV_markers, file = "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/DEGs/BAP1vsV_all_markers.csv")

PBRM1vV_markers <- FindMarkers(PDOT, ident.1 = "PBRM1_Driven", ident.2 = "VHL_monodriver", max.cells.per.ident = 500, min.pct = 0.25, logfc.threshold = 0.5)
PBRM1vV_markers <- PBRM1vV_markers[PBRM1vV_markers$p_val_adj < 0.05,]
write.csv(PBRM1vV_markers, file = "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/DEGs/PBRM1vsV_all_markers.csv")

# markers BAP1 vs. PBRM1
BAP1vsPBRM1_markers <- FindMarkers(PDOT, ident.1 = "BAP1_Driven", ident.2 = "PBRM1_Driven", min.pct = 0.25, logfc.threshold = 0.5)
BAP1vsPBRM1_markers <- BAP1vsPBRM1_markers[BAP1vsPBRM1_markers$p_val_adj < 0.05,]
write.csv(BAP1vsPBRM1_markers, file = "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/DEGs/BAP1vsPBRM1_markers.csv")


## other comparisons - co-occurence compared to normal or monodriver
BAP1_PBRM1vN_markers <- FindMarkers(PDOT, ident.1 = "PBRM1_BAP1", ident.2 = "Normal", max.cells.per.ident = 500, min.pct = 0.25, logfc.threshold = 0.5)
BAP1_PBRM1vN_markers <- BAP1_PBRM1vN_markers[BAP1_PBRM1vN_markers$p_val_adj < 0.05,]
write.csv(BAP1_PBRM1vN_markers, file = "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/DEGs/BAP1_PBRM1vsN_markers.csv")

BAP1_PBRM1vV_markers <- FindMarkers(PDOT, ident.1 = "PBRM1_BAP1", ident.2 = "VHL_monodriver", max.cells.per.ident = 500, min.pct = 0.25, logfc.threshold = 0.5)
BAP1_PBRM1vV_markers <- BAP1_PBRM1vV_markers[BAP1_PBRM1vV_markers$p_val_adj < 0.05,]
write.csv(BAP1_PBRM1vV_markers, file = "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/DEGs/BAP1_PBRM1vsV_markers.csv")

[1] "K1013_GB1" "K1022_GB1" "K1030_GB1" "K1030_GB2" "K1030_GN1" "K1040_GB1" "K1042_GB1" "K1042_GE2" "K1063_GB1"
[10] "K1063_GB2" "K1184_GA1" "K1184_GC1"
## other comparisons - find markers between normal and tumour groups
PDOT <- SetIdent(PDOT, value = PDOT@meta.data$LineID)
# Define pairs of comparisons, compared to N
ident_pairs <- list(
  "K1013_GB1" = "K1030_GN1",
  "K1022_GB1" = "K1030_GN1",
  "K1030_GB1" = "K1030_GN1",
  "K1030_GB2" = "K1030_GN1",
  "K1040_GB1" = "K1030_GN1",
  "K1042_GB1" = "K1030_GN1",
  "K1042_GE2" = "K1030_GN1",
  "K1063_GB1" = "K1030_GN1",
  "K1063_GB2" = "K1030_GN1",
  "K1184_GA1" = "K1030_GN1",
  "K1184_GC1" = "K1030_GN1"
)
# Define the output directory
output_dir <- "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/DEGs/"
# Loop through ident_pairs
for (ident.1 in names(ident_pairs)) {
  ident.2 <- ident_pairs[[ident.1]]
  
  # Run FindMarkers
  markers <- FindMarkers(
    PDOT, 
    ident.1 = ident.1, 
    ident.2 = ident.2, 
    max.cells.per.ident = 1200, 
    min.pct = 0.25, 
    logfc.threshold = 0.5
  )
  
  # Filter by adjusted p-value
  markers_filtered <- markers[markers$p_val_adj < 0.05, ]
  # Create a filename
  output_file <- paste0(output_dir, ident.1, "vs", ident.2, ".csv")
  
  # Save to CSV
  write.csv(markers_filtered, file = output_file, row.names = TRUE)
}
# Define pairs of comparisons, compared to V
ident_pairs <- list(
  "K1022_GB1" = "K1013_GB1",
  "K1030_GB1" = "K1013_GB1",
  "K1030_GB2" = "K1013_GB1",
  "K1040_GB1" = "K1013_GB1",
  "K1042_GB1" = "K1013_GB1",
  "K1042_GE2" = "K1013_GB1",
  "K1063_GB1" = "K1013_GB1",
  "K1063_GB2" = "K1013_GB1",
  "K1184_GA1" = "K1013_GB1",
  "K1184_GC1" = "K1013_GB1"
)
# Define the output directory
output_dir <- "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/DEGs/"
# Loop through ident_pairs
for (ident.1 in names(ident_pairs)) {
  ident.2 <- ident_pairs[[ident.1]]
  
  # Run FindMarkers
  markers <- FindMarkers(
    PDOT, 
    ident.1 = ident.1, 
    ident.2 = ident.2, 
    max.cells.per.ident = 1200, 
    min.pct = 0.25, 
    logfc.threshold = 0.5
  )
  
  # Filter by adjusted p-value
  markers_filtered <- markers[markers$p_val_adj < 0.05, ]
  # Create a filename
  output_file <- paste0(output_dir, ident.1, "vs", ident.2, ".csv")
  
  # Save to CSV
  write.csv(markers_filtered, file = output_file, row.names = TRUE)
}
########################################################################## 
# Pathway Enrichment Analysis - done in local 
##########################################################################
PDOT <- readRDS("/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/Seurats/PDOT_v2_clean_int.rds")


########################################################################## 
# NMF deconvolution
##########################################################################
PDOT <- readRDS("/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/Seurats/PDOT_v3_correct.rds")
# PDOT <- readRDS("/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/Seurats/PDOT_v2_clean_int.rds")

DefaultAssay(PDOT)
seu.list <- SplitObject(PDOT, split.by = "LineID")

geneNMF.programs <- multiNMF(seu.list, assay="RNA", k=4:9, min.exp = 0.05)

geneNMF.metaprograms <- getMetaPrograms(geneNMF.programs,
                                        metric = "cosine",
                                        weight.explained = 0.7,
                                        nMP=10,
                                        min.confidence = 0.6)

ph <- plotMetaPrograms(geneNMF.metaprograms,
                       similarity.cutoff = c(0.1,1))
ph

geneNMF.metaprograms$metaprograms.metrics
lapply(geneNMF.metaprograms$metaprograms.genes, head)

mp.genes <- geneNMF.metaprograms$metaprograms.genes
seu <- AddModuleScore_UCell(PDOT, features = mp.genes, ncores=4, name = "")
VlnPlot(seu, features=names(mp.genes), group.by = "Evo_Type",
        pt.size = 0, ncol=5)
top_p <- lapply(geneNMF.metaprograms$metaprograms.genes, function(program) {
  runGSEA(program, universe=rownames(PDOT), category = "C5", subcategory = "GO:BP")
})

########################################################################## 
# WGCNA deconvolution
##########################################################################



########################################################################## 
# SCENIC TF activity inference - done in python and bash, read in meta
##########################################################################
library(SCENIC)
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCopeLoomR)
# library(loomR)
library(SeuratDisk)
setwd("./SCENIC/")
### get expression matrix and save loom
PDOT <- readRDS("/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/Seurats/PDOT_v4_gms.rds")
exprMat <- GetAssayData(PDOT, layer = "counts")
cellInfo <- data.frame(seuratCluster=Idents(PDOT))
cellInfo <- data.frame(cellInfo)
cellInfo$LineID <- PDOT@meta.data[["LineID"]]
cellInfo$Evo_Type <- PDOT@meta.data[["Evo_Type"]]
loom <- build_loom("./data/PDOT.loom", dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)

### Load data
loom <- open_loom("./data/PDOT.loom")
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)

### SAVE Cell info/metadata
dir.create("./int")
saveRDS(cellInfo, file="./int/cellInfo.Rds")

### LOOK AT the results
PDOT <- readRDS("/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/Seurats/PDOT_v4_gms.rds")
# read the dimension reduction from SCENIC auc_matrix and add to PDOT
scenic_umap <- read.delim("/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/SCENIC/output/RUN_1/scenic_umap.txt", row.names=1)
colnames(scenic_umap) <- c("scenicumap_1", "scenicumap_2")
all(rownames(scenic_umap) == Cells(PDOT))
PDOT[["scenic_umap"]] <- CreateDimReducObject(embeddings = as.matrix(scenic_umap), key = "scenicumap_", assay = DefaultAssay(PDOT))
DimPlot(PDOT, reduction = "scenic_umap")
DP_Evo_Type_regulon <- DimPlot(PDOT, reduction = "scenic_umap", 
                       group.by = "Evo_Type", 
                       label = T, 
                       label.size = 5, 
                       pt.size = 1.5, 
                       repel = TRUE,
                       shuffle = T) +
  ggplot2::ggtitle("UMAP Projection") +
  ggplot2::theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text = element_blank(),      # Remove axis text
    axis.title = element_blank(),     # Remove axis titles
    axis.ticks = element_blank(),     # Remove axis ticks
    legend.position = "right",
    legend.text = element_text(size = 10)
  )
ggsave(filename = "/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/Results/Plots/DimPlot/DP_Evo_Type_regulon.pdf", plot = DP_Evo_Type_regulon, width = 8, height = 6, dpi = 300)

auc_mtx_1 <- read.csv("/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/SCENIC/output/RUN_1/auc_mtx.csv", row.names=1)
colnames(auc_mtx_1) <- substr(colnames(auc_mtx_1), 1, nchar(colnames(auc_mtx_1)) - 3)
# Calculate average for each group
auc_mtx_1$Evo_Type <- cellInfo$Evo_Type
group_means <- auc_mtx_1 %>%
  group_by(Evo_Type) %>%
  summarise(across(everything(), mean))
auc_mtx_zscore_1 <- as.data.frame(group_means[-1])
rownames(auc_mtx_zscore_1) <- group_means$Evo_Type
auc_mtx_zscore_1 <- scale(auc_mtx_zscore_1)
auc_mtx_zscore_1 <- auc_mtx_zscore_1[c("Normal", "VHL_monodriver", "BAP1_Driven", "PBRM1_Driven", "BAP1_PBRM1"), ]


auc_mtx_2 <- read.csv("/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/SCENIC/output/RUN_2/auc_mtx.csv", row.names=1)
colnames(auc_mtx_2) <- substr(colnames(auc_mtx_2), 1, nchar(colnames(auc_mtx_2)) - 3)
# Calculate average for each group
auc_mtx_2$Evo_Type <- cellInfo$Evo_Type
group_means <- auc_mtx_2 %>%
  group_by(Evo_Type) %>%
  summarise(across(everything(), mean))
auc_mtx_zscore_2 <- as.data.frame(group_means[-1])
rownames(auc_mtx_zscore_2) <- group_means$Evo_Type
auc_mtx_zscore_2 <- scale(auc_mtx_zscore_2)
auc_mtx_zscore_2 <- auc_mtx_zscore_2[c("Normal", "VHL_monodriver", "BAP1_Driven", "PBRM1_Driven", "BAP1_PBRM1"), ]


auc_mtx_3 <- read.csv("/nemo/project/proj-tracerX/working/RENALSC/DAQI/Scratch/PDOT_Parse/SCENIC/output/RUN_3/auc_mtx.csv", row.names=1)
colnames(auc_mtx_3) <- substr(colnames(auc_mtx_3), 1, nchar(colnames(auc_mtx_3)) - 3)
# Calculate average for each group
auc_mtx_3$Evo_Type <- cellInfo$Evo_Type
group_means <- auc_mtx_3 %>%
  group_by(Evo_Type) %>%
  summarise(across(everything(), mean))
auc_mtx_zscore_3 <- as.data.frame(group_means[-1])
rownames(auc_mtx_zscore_3) <- group_means$Evo_Type
auc_mtx_zscore_3 <- scale(auc_mtx_zscore_3)
auc_mtx_zscore_3 <- auc_mtx_zscore_3[c("Normal", "VHL_monodriver", "BAP1_Driven", "PBRM1_Driven", "BAP1_PBRM1"), ]


overlap <- intersect(colnames(auc_mtx_zscore_1), intersect(colnames(auc_mtx_zscore_2), colnames(auc_mtx_zscore_3)))
auc_mtx_zscore_1 <- auc_mtx_zscore_1[, overlap]
auc_mtx_zscore_2 <- auc_mtx_zscore_2[, overlap]
auc_mtx_zscore_3 <- auc_mtx_zscore_3[, overlap]
auc_mtx_zscore <- auc_mtx_zscore_1[,top]
pheatmap(auc_mtx_zscore, cluster_cols = T, cluster_rows = F)

consistent_list_1 <- c("AR", "ATF3", "BACH1", "BCLAF1", "BPTF", "BRCA1", "CREB5", "E2F1", "E2F2", "E2F3", "E2F7", "E2F8", "EBF1", "ELF1", "ELF2", "ELF5", "ELK4", "EP300", "ETS2", "FOS", "FOSB", "FOXE1", "FOXF1", "GATA3", "HIVEP1", "HMGA2", "HOXA13", "IKZF2", "JUN", "JUNB", "KLF6", "MBD2", "MYB", "MYBL1", "NFATC1", "NFE2L1", "NFE2L2", "NFIA", "NFIB", "NFKB1", "NFKB2", "NR2F1", "NR3C1", "PBX1", "RARG", "RB1", "RFX2", "RFX3", "RUNX3", "SP3", "SREBF2", "STAT1", "STAT2", "STAT6", "TCF4", "TCF7L1", "TFDP1", "THRB", "ZEB1", "YBX1", "ZNF354C")

consistent_list <- c("AR", "ATF3", "ATF4", "BACH1", "BACH2", "BCLAF1", "BHLHE40", "BPTF", "BRCA1", "CEBPA", "CREB5", "E2F1", "E2F2", "E2F3", "E2F7", "E2F8", "EBF1", "ELF1", "ELF5", "ELK4", "EP300", "ETS2", "ETV1", "FLI1", "FOS", "FOSB", "FOXA1", "FOXE1", "FOXF1", "GATA3", "HIVEP1", "HMGA2", "HNF1A", "HNF4A", "HNF4G", "HOXA13", "HOXB2", "HOXB6", "HOXD10", "IKZF2", "IRF1", "IRF7", "JUN", "JUNB", "JUND", "KLF6", "LHX9", "MBD2", "MYB", "MYBL1", "NFATC1", "NFE2L1", "NFE2L2", "NFIA", "NFIB", "NFKB1", "NFKB2", "NR2F1", "NR2F2", "NR3C1", "PAX6", "PBX1", "POU3F1", "POU5F1", "RB1", "REL", "RFX2", "RFX3", "RUNX3", "SMARCC2", "SP3", "SPDEF", "SREBF2", "SRF", "STAT1", "STAT2", "STAT3", "STAT6", "TCF4", "TCF7L1", "TFDP1", "UBE2K", "ZEB1", "ZNF148", "ZNF354C")
TF_list <- c("IKZF2", "ELK1", "ELF5", "GATA3", "REL", "ZEB1", "KLF6", "ETS1", "IRF1", "STAT6", "NFKB1", "NFKB2", "FOS", "FOSB", "JUN", "JUND", "JUNB", "ELK4", "E2F1", "BRCA1", "BHLHE40", "TFE3", "HMGA2", "ETV4", "HOXA9", "NFIB", "NFIA", "HDX", "TEAD4", "THRB", "BACH1", "FOXA1", "GATA6", "FLI1", "ZNF354C", "ZNF528", "MYB", "TCF4", "TCF7L1")
auc_mtx_zscore <- auc_mtx_zscore_3[, overlap]


auc_mtx_1$Evo_Type <- cellInfo$LineID
group_means <- auc_mtx_1 %>%
  group_by(Evo_Type) %>%
  summarise(across(everything(), mean))
auc_mtx_zscore_1 <- as.data.frame(group_means[-1])
rownames(auc_mtx_zscore_1) <- group_means$Evo_Type
auc_mtx_zscore_1 <- scale(auc_mtx_zscore_1)
write.csv(auc_mtx_zscore_1, file = "./SCENIC/output/RUN_1/Regulon_zscore.csv")

auc_mtx_zscore_1 <- auc_mtx_zscore_1[c("Normal", "VHL_monodriver", "BAP1_Driven", "PBRM1_Driven", "BAP1_PBRM1"), ]
auc_mtx_zscore_1 <- auc_mtx_zscore_1[c("K1030_GN1", "K1013_GB1", "K1042_GE2", "K1042_GB1", "K1184_GA1", "K1184_GC1", "K1040_GB1", "K1030_GB2", "K1030_GB1", "K1063_GB1", "K1063_GB2", "K1022_GB1"),]


# prioritse regulon
for (x in overlap) {
  a <- pheatmap(as.data.frame(cbind(auc_mtx_zscore_1[,x], auc_mtx_zscore_2[,x], auc_mtx_zscore_3[,x])))
  ggsave(filename = paste0("./temp/hm_", x, ".pdf"), a)
}

vars <- apply(auc_mtx_1, 2, var, na.rm = TRUE)
# Sort descending and take top 20
top100_cols <- names(sort(vars, decreasing = TRUE))[1:100]
# View them
print(top100_cols)
top <- intersect(top100_cols, colnames(auc_mtx_zscore_1))
top <- top[1:60]
pheatmap(auc_mtx_zscore_1, cluster_cols = T, cluster_rows = F)
