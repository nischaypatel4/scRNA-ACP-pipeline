#---------------------------------------------
#---------------------------------------------
# ACP Scimilarity Input Manipulation
# Input : acp2_final
# Output: acp2_epi_cells_raw_counts.csv and acp2_epi_cells_metadata.csv (Input to Scimilarity Tool)
#---------------------------------------------
#---------------------------------------------

#---------------------------------------------
# 0. Load libraries
#---------------------------------------------
library(Seurat)         # Core single-cell toolkit
library(ggplot2)        # Data visualization
library(patchwork)      # Combine plots
library(sctransform)    # SCTransform normalization
library(harmony)        # Batch correction (not used later but loaded)
library(cowplot)        # Plotting utilities

#---------------------------------------------
# Set directories for processed data and Scimilarity input files
#---------------------------------------------
PROC_DIR <- "data/processed"
OUT_DIR  <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/Scimilarity_Input"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)   # Create output folder if missing

#---------------------------------------------
# 1. Load final filtered objects
#---------------------------------------------
acp2 <- readRDS(file.path(PROC_DIR, "acp2_final.rds"))   # Load Seurat object for ACP2

#---------------------------------------------
# ACP2: Initial Preprocessing and Clustering
#---------------------------------------------

# Normalize raw counts (default log-normalization)
acp2_swf <- NormalizeData(acp2)

# Identify highly variable genes
acp2_swf <- FindVariableFeatures(acp2_swf)

# Scale data (centering and variance scaling)
acp2_swf <- ScaleData(acp2_swf)

# Perform PCA on scaled data
acp2_swf <- RunPCA(acp2_swf, ndims = 50, verbose = FALSE)

# Visualize variance explained to pick dimensions
Elbow_acp2 <- ElbowPlot(acp2_swf, ndims = 50, reduction = "pca")
Elbow_acp2

# Build nearest neighbor graph and cluster cells at multiple resolutions
acp2_swf <- FindNeighbors(acp2_swf, dims = 1:25, reduction = "pca")
acp2_swf <- FindClusters(acp2_swf, resolution = c(0.1, 0.2, 0.3, 0.5, 0.7, 1))

# Generate UMAP for visual exploration
acp2_swf <- RunUMAP(acp2_swf, dims = 1:25, reduction = "pca")

# Open metadata in viewer (interactive RStudio function)
View(acp2_swf@meta.data)

# Check identities of clusters and assign clustering resolution
Idents(acp2_swf)
Idents(acp2_swf) <- "RNA_snn_res.0.2"

# Plot UMAP by cluster assignments
Plot_acp2 <- DimPlot(acp2_swf, group.by = c("RNA_snn_res.0.2"), ncol = 1, reduction = "umap")
Plot_acp2

#---------------------------------------------
# Marker Gene Expression Visualization
#---------------------------------------------

# Create UMAP FeaturePlots for major lineage marker genes
p1_acp2 <- FeaturePlot(acp2_swf, features = c("CDH1", "EPCAM"), reduction = "umap", combine = FALSE)
p2_acp2 <- FeaturePlot(acp2_swf, features = c("LYZ", "CD68"), reduction = "umap", combine = FALSE)
p3_acp2 <- FeaturePlot(acp2_swf, features = c("CD3E"), reduction = "umap", combine = FALSE)
p4_acp2 <- FeaturePlot(acp2_swf, features = c("CD79A"), reduction = "umap", combine = FALSE)
p5_acp2 <- FeaturePlot(acp2_swf, features = c("AQP4", "GFAP"), reduction = "umap", combine = FALSE)
p6_acp2 <- FeaturePlot(acp2_swf, features = c("SOX10", "OLIG2"), reduction = "umap", combine = FALSE)
p7_acp2 <- FeaturePlot(acp2_swf, features = c("MOG"), reduction = "umap", combine = FALSE)
p8_acp2 <- FeaturePlot(acp2_swf, features = c("PECAM1"), reduction = "umap", combine = FALSE)

# Combine all marker plots into a single multi-panel plot
all_plots_acp2 <- c(p1_acp2, p2_acp2, p3_acp2, p4_acp2, p5_acp2, p6_acp2, p7_acp2, p8_acp2)
combined_plot_acp2 <- plot_grid(plotlist = all_plots_acp2, ncol = 3)

# Display marker expression and cluster plot side by side
print(combined_plot_acp2)
combined_plot_acp2 | Plot_acp2

#---------------------------------------------
# Annotate Cell Types from Clusters
#---------------------------------------------

# Create mapping from clusters to known cell types
cluster_ids_acp2 <- as.character(acp2_swf$RNA_snn_res.0.2)

celltype_annotation_acp2 <- dplyr::case_when(
  cluster_ids_acp2 %in% c("1", "3", "4") ~ "Epithelial cells",
  cluster_ids_acp2 == "0" ~ "Myeloid cells",
  cluster_ids_acp2 == "6" ~ "T cells",
  cluster_ids_acp2 == "8" ~ "B cells",
  cluster_ids_acp2 == "5" ~ "Astrocyte",
  cluster_ids_acp2 == "2" ~ "OPC",
  cluster_ids_acp2 == "7" ~ "Oligodendrocyte",
  
  TRUE ~ "Unknown"
)

# Store cell type annotations in metadata
acp2_swf$celltype <- celltype_annotation_acp2

# Plot annotated UMAP with colors for each cell type
Annotation_acp2 <- DimPlot(
  acp2_swf,
  group.by = "celltype",
  reduction = "umap",
  label = FALSE,
  cols = c(
    "Epithelial cells" = "blue",
    "Myeloid cells" = "orange",
    "T cells" = "red",
    "B cells" = "yellow",
    "Astrocyte" = "pink",
    "OPC" = "black",
    "Oligodendrocyte" = "green",
    "Unknown" = "grey"
  )
)
Annotation_acp2 
Annotation_acp2 | combined_plot_acp2

#---------------------------------------------
# Subset and Re-cluster Epithelial Cells
#---------------------------------------------

# Extract only epithelial cells for deeper analysis
epi_cells_acp2 <- subset(acp2_swf, subset = celltype == "Epithelial cells")
View(epi_cells_acp2@meta.data)

# Run standard Seurat workflow on epithelial subset
epi_cells_acp2 <- FindVariableFeatures(epi_cells_acp2)
epi_cells_acp2 <- ScaleData(epi_cells_acp2)
epi_cells_acp2 <- RunPCA(epi_cells_acp2, ndims = 50, verbose = FALSE)
ElbowPlot(epi_cells_acp2, ndims = 50, reduction = "pca")
epi_cells_acp2 <- FindNeighbors(epi_cells_acp2, dims = 1:25, reduction = "pca")
epi_cells_acp2 <- FindClusters(epi_cells_acp2, resolution = c(0.1, 0.2, 0.3, 0.5, 0.7, 1))
epi_cells_acp2 <- RunUMAP(epi_cells_acp2, dims = 1:25, reduction = "pca")

# Set cluster identities for epithelial cells
Idents(epi_cells_acp2) <- "RNA_snn_res.0.1"
Plot_epi_cells_acp2 <- DimPlot(epi_cells_acp2, group.by = "RNA_snn_res.0.1", reduction = "umap", ncol = 1)
Plot_epi_cells_acp2 | Annotation_acp2

#---------------------------------------------
# Module Scores: Define Gene Signatures
#---------------------------------------------

# Define four epithelial signatures (CC, WE, KE, PE)
CC_genes <- c("KIAA0101", "TOP2A", "CDK1", "UBE2C", "NUSAP1", "PBK", "CENPM", "TPX2",
              "BIRC5", "ZWINT", "CCNA2", "ESCO2", "CENPF", "NCAPH", "MKI67", "CDCA3",
              "TK1", "NDC80", "CKAP2L", "KIFC1", "NCAPG", "SGOL1", "HMMR", "NUF2",
              "AURKB", "CEP55", "CENPK", "SPC25", "CDKN3", "CDCA5")

WE_genes <- c("MUCL1", "RBMS3", "COL1A2", "KRT23", "VCAN", "FGF4", "DEFA5", "TANGO6",
              "FGF9", "FGF20", "DPEP1", "SEZ6L", "WNT5A", "WNT10A", "NKX2-5", "SLC1A5",
              "SEMA3C", "CAMK2B", "EFHD1", "KREMEN2", "LMO2", "NOTUM", "WNT10B", "TMEM45B",
              "DKK4", "RGPD3", "ITGA1", "PCCA", "CTNNB1", "MUC1")

KE_genes <- c("SCUBE3", "TAGLN", "KIF5C", "CALB1", "ISG15", "CD24", "SYTL2", "CALML3",
              "ODAM", "HES1", "ATP6V0D2", "KRT75", "CLDN10", "MGP", "IFI6", "TM4SF1",
              "MACC1", "KRT6B", "VIM", "HHIP", "TRIB1", "PROCR", "GBP1", "MX1",
              "KRT16", "GJB6", "SLC20A2", "KLF6", "KRT6C", "IFIT1", "MAF", "FLNB",
              "GJB2", "CLDN1", "MYC", "KRT7")

PE_genes <- c("CTNNAL1", "SERTAD1", "CRABP1", "MT1E", "ATF3", "DLX6-AS1", "NFKBIA", "MT1X",
              "CYR61", "IFRD1", "EGR1", "MT1F", "NR4A1", "JUNB", "SLC39A10", "TSPAN8",
              "FOS", "MT1G", "KLF10", "CD47", "SFRP1", "FJX1", "EGR3", "EPHA7",
              "IRF1", "WIF1", "COL17A1", "FOSB", "FRZB", "IER3", "JUN", "JUND", "HSPA1A",
              "EGR2", "DDIT3", "POSTN", "ARL4D", "FAM41C", "SLC35G1", "KCTD12", "EPB41L2",
              "GADD45B", "PTCH1", "KLHL42", "VSTM2B")

#---------------------------------------------
# Calculate Module Scores for Each Cell
#---------------------------------------------

# Combine all layers to make sure module score calculation works
epi_cells_acp2_1 <- JoinLayers(epi_cells_acp2, overwrite = TRUE)

# Compute module scores separately for each signature
epi_cells_acp2_1 <- AddModuleScore(epi_cells_acp2_1, features = list(CC_genes), name = "CC_Score")
epi_cells_acp2_1 <- AddModuleScore(epi_cells_acp2_1, features = list(WE_genes), name = "WE_Score")
epi_cells_acp2_1 <- AddModuleScore(epi_cells_acp2_1, features = list(KE_genes), name = "KE_Score")
epi_cells_acp2_1 <- AddModuleScore(epi_cells_acp2_1, features = list(PE_genes), name = "PE_Score")

View(epi_cells_acp2@meta.data)
View(epi_cells_acp2_1@meta.data)

# Plot module scores on UMAP
epi_featureplot_acp2 <- FeaturePlot(
  epi_cells_acp2_1,
  features = c("CC_Score1", "WE_Score1", "KE_Score1", "PE_Score1"),
  ncol = 2
) + plot_annotation(title = "Acp2 Module Scores on Epithelial UMAP")
epi_featureplot_acp2

#---------------------------------------------
# Assign Dominant Epithelial Signature per Cell
#---------------------------------------------

# Determine which signature scores highest per cell
score_cols <- c("CC_Score1", "WE_Score1", "KE_Score1", "PE_Score1")
scores <- epi_cells_acp2_1@meta.data[, score_cols]
best_idx_acp2 <- apply(scores, 1, which.max)
epi_cells_acp2_1$epi_type <- score_cols[best_idx_acp2]

# Set colors for module identity visualization
epi_colors <- c(
  CC_Score1 = "lightblue",
  WE_Score1 = "darkgreen",
  KE_Score1 = "orange",
  PE_Score1 = "purple"
)

# Plot dominant epithelial signature identity
epi_umap_acp2 <- DimPlot(
  epi_cells_acp2_1,
  group.by = "epi_type",
  reduction = "umap",
  label = FALSE
) +
  ggtitle("Dominant Epithelial Signature per Cell") +
  scale_color_manual(values = epi_colors)

epi_umap_acp2

View(epi_cells_acp2_1@meta.data)

#---------------------------------------------
# Annotate Subclusters within Epithelial Cells
#---------------------------------------------

# Extract cluster IDs for epithelial cells
cluster_ids_epi_acp2 <- as.character(epi_cells_acp2$RNA_snn_res.0.1)

# Map subclusters to epithelial types
celltype_annotation_acp2 <- dplyr::case_when(
  cluster_ids_epi_acp2 == "1" ~ "PE",
  cluster_ids_epi_acp2 == "2" ~ "KE",
  cluster_ids_epi_acp2 == "4" ~ "WE",
)

# Assign annotations back to epithelial subset
epi_cells_acp2$celltype <- celltype_annotation_acp2

View(epi_cells_acp2@meta.data)

# Plot UMAP with epithelial subtypes labeled
Final_acp2 <- DimPlot(
  epi_cells_acp2,
  group.by = "celltype",
  reduction = "umap",
  label = FALSE,
  cols = c(
    "PE" = "blue",
    "KE" = "orange",
    "WE" = "darkgreen",
    "CC" = "pink",
    "Unknown" = "grey"
  )
)
Final_acp2

#---------------------------------------------
# Save Raw Counts and Metadata
#---------------------------------------------

# Pull out raw counts (non-log normalized)
counts_acp2 <- GetAssayData(epi_cells_acp2, assay = "RNA", slot = "counts")
counts_df_acp2 <- as.data.frame(as.matrix(counts_acp2))

View(counts_acp2)
View(counts_df_acp2)

# Extract metadata and attach cell_id column
metadata_df_acp2 <- epi_cells_acp2@meta.data
metadata_df_acp2$cell_id <- rownames(metadata_df_acp2)

# Save counts and metadata for Scimilarity input
write.csv(counts_df_acp2, file = file.path(OUT_DIR, "acp2_epi_cells_raw_counts.csv"))
write.csv(metadata_df_acp2, file = file.path(OUT_DIR, "acp2_epi_cells_metadata.csv"), row.names = FALSE)

# Plot combination of clustering and module-based annotation
Plot_epi_cells_acp2 + epi_umap_acp2
Final_acp2 + epi_umap_acp2