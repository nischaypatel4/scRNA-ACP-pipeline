#---------------------------------------------
#---------------------------------------------
# ACP Scimilarity Input Manipulation
# Input : acp1_final/acp2_final/acp3_final
# Output: acp1/2/3_epi_cells_raw_counts.csv and acp1/2/3_epi_cells_metadata.csv (Input to Scimilarity Tool)
#---------------------------------------------
#---------------------------------------------

#---------------------------------------------
# 0. Load libraries
#---------------------------------------------
library(Seurat)          # Load Seurat package for single-cell RNA-seq data analysis
library(ggplot2)         # Load ggplot2 for plotting
library(patchwork)       # Load patchwork to combine ggplot2 plots
library(sctransform)     # Load sctransform package (normalization, variance stabilization)
library(harmony)         # Load harmony for batch correction
library(cowplot)         # Load cowplot for combining multiple plots

#---------------------------------------------
# Set directories
#---------------------------------------------
PROC_DIR <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed"  # Directory with processed data
OUT_DIR  <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/Scimilarity_Input"  # Output directory for Scimilarity input
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)  # Create output directory if it doesn't exist, suppress warnings

#---------------------------------------------
# 1. Load final filtered objects
#---------------------------------------------
acp1 <- readRDS(file.path(PROC_DIR, "acp1_final.rds"))  # Load the processed Seurat object for acp1

#---------------------------------------------
# ACP1: Initial Preprocessing and Clustering
#---------------------------------------------

# Normalize raw counts using default log-normalization
acp1_swf <- NormalizeData(acp1)

# Identify variable features (highly variable genes)
acp1_swf <- FindVariableFeatures(acp1_swf)

# Scale data (center and scale expression values)
acp1_swf <- ScaleData(acp1_swf)

# Perform Principal Component Analysis (PCA) with 50 components
acp1_swf <- RunPCA(acp1_swf, ndims = 50, verbose = FALSE)

# Construct neighbor graph using first 25 PCs and cluster cells at multiple resolutions
acp1_swf <- FindNeighbors(acp1_swf, dims = 1:25, reduction = "pca")
acp1_swf <- FindClusters(acp1_swf, resolution = c(0.1, 0.2, 0.3, 0.5, 0.7, 1))

# Run UMAP dimensionality reduction for visualization
acp1_swf <- RunUMAP(acp1_swf, dims = 1:25, reduction = "pca")

# Check current identities of cells
Idents(acp1_swf)

# Set cluster identities to the resolution 0.2 clustering
Idents(acp1_swf) <- "RNA_snn_res.0.2"

# Plot UMAP with clusters colored by resolution 0.2
Plot_acp1 <- DimPlot(acp1_swf, group.by = c("RNA_snn_res.0.2"), ncol = 1, reduction = "umap")

#---------------------------------------------
# Marker Gene Expression Visualization
#---------------------------------------------
# Generate UMAP feature plots for known marker genes (markers for cell types)

p1_acp1 <- FeaturePlot(acp1_swf, features = c("CDH1", "EPCAM"), reduction = "umap", combine = FALSE)  # Epithelial markers
p2_acp1 <- FeaturePlot(acp1_swf, features = c("LYZ", "CD68"), reduction = "umap", combine = FALSE)    # Myeloid markers
p3_acp1 <- FeaturePlot(acp1_swf, features = c("CD3E"), reduction = "umap", combine = FALSE)           # T cell marker
p4_acp1 <- FeaturePlot(acp1_swf, features = c("CD79A"), reduction = "umap", combine = FALSE)          # B cell marker
p5_acp1 <- FeaturePlot(acp1_swf, features = c("AQP4", "GFAP"), reduction = "umap", combine = FALSE)   # Astrocyte markers
# p6_acp1 <- FeaturePlot(acp1_swf, features = c("SOX10", "OLIG2", "GFAP"), reduction = "umap", combine = FALSE)  # Commented out: oligodendrocyte markers
p7_acp1 <- FeaturePlot(acp1_swf, features = c("PECAM1"), reduction = "umap", combine = FALSE)         # Endothelial marker

# Combine all marker plots into a grid layout with 3 columns
all_plots_acp1 <- c(p1_acp1, p2_acp1, p3_acp1, p4_acp1, p5_acp1, p7_acp1)
combined_plot_acp1 <- plot_grid(plotlist = all_plots_acp1, ncol = 3)

# Display combined marker gene expression plots
print(combined_plot_acp1)

# Display combined marker plots side-by-side with cluster UMAP
combined_plot_acp1 | Plot_acp1

#---------------------------------------------
# Annotate Cell Types from Clusters
#---------------------------------------------

# Extract cluster IDs as character vector from clustering at resolution 0.2
cluster_ids_acp1 <- as.character(acp1_swf$RNA_snn_res.0.2)

# Assign known cell type labels to clusters using case_when
celltype_annotation_acp1 <- dplyr::case_when(
  cluster_ids_acp1 %in% c("0", "1", "2", "3", "5", "6") ~ "Epithelial cells",  # clusters for epithelial cells
  cluster_ids_acp1 == "4" ~ "Myeloid cells",                                  # cluster 4 as myeloid cells
  cluster_ids_acp1 == "8" ~ "T cells",                                        # cluster 8 as T cells
  cluster_ids_acp1 == "7" ~ "Astrocyte",                                      # cluster 7 as astrocytes
  TRUE ~ "Unknown"                                                            # everything else unknown
)

# Add cell type annotation as metadata column
acp1_swf$celltype <- celltype_annotation_acp1

# Plot UMAP with cells colored by annotated cell type
Annotation_acp1 <- DimPlot(
  acp1_swf,
  group.by = "celltype",   # group by cell type annotation
  reduction = "umap",      # use UMAP embedding
  label = FALSE,           # do not add cluster labels
  cols = c(                # assign colors for cell types
    "Epithelial cells" = "blue",
    "Myeloid cells" = "orange",
    "T cells" = "red",
    "Astrocyte" = "pink",
    "Unknown" = "grey"
  )
)

# Display annotation plot
Annotation_acp1 

# Show side-by-side with combined marker plots
Annotation_acp1 | combined_plot_acp1

#---------------------------------------------
# Subset and Re-cluster Epithelial Cells
#---------------------------------------------

# Subset only epithelial cells from the Seurat object
epi_cells_acp1 <- subset(acp1_swf, subset = celltype == "Epithelial cells")

# Identify variable features in epithelial subset
epi_cells_acp1 <- FindVariableFeatures(epi_cells_acp1)

# Scale data in epithelial subset
epi_cells_acp1 <- ScaleData(epi_cells_acp1)

# Run PCA with 50 components on epithelial subset
epi_cells_acp1 <- RunPCA(epi_cells_acp1, ndims = 50, verbose = FALSE)

# Plot elbow plot to decide how many PCs to use
ElbowPlot(epi_cells_acp1, ndims = 50, reduction = "pca")

# Find neighbors using PCs 1-25
epi_cells_acp1 <- FindNeighbors(epi_cells_acp1, dims = 1:25, reduction = "pca")

# Cluster epithelial cells at multiple resolutions
epi_cells_acp1 <- FindClusters(epi_cells_acp1, resolution = c(0.1, 0.2, 0.3, 0.5, 0.7, 1))

# Run UMAP on epithelial cells using PCs 1-25
epi_cells_acp1 <- RunUMAP(epi_cells_acp1, dims = 1:25, reduction = "pca")

# Set cluster identity to resolution 0.1 clustering for epithelial cells
Idents(epi_cells_acp1) <- "RNA_snn_res.0.1"

# Plot UMAP of epithelial cells colored by clusters at resolution 0.1
Plot_epi_cells_acp1 <- DimPlot(epi_cells_acp1, group.by = "RNA_snn_res.0.1", reduction = "umap", ncol = 1)

# Display epithelial clusters alongside previous annotation UMAP
Plot_epi_cells_acp1 | Annotation_acp1

#---------------------------------------------
# Module Scores: Define Gene Signatures
#---------------------------------------------

# Define gene sets corresponding to epithelial subtypes

CC_genes <- c("KIAA0101", "TOP2A", "CDK1", "UBE2C", "NUSAP1", "PBK", "CENPM", "TPX2",
              "BIRC5", "ZWINT", "CCNA2", "ESCO2", "CENPF", "NCAPH", "MKI67", "CDCA3",
              "TK1", "NDC80", "CKAP2L", "KIFC1", "NCAPG", "SGOL1", "HMMR", "NUF2",
              "AURKB", "CEP55", "CENPK", "SPC25", "CDKN3", "CDCA5")  # Cell cycle genes

WE_genes <- c("MUCL1", "RBMS3", "COL1A2", "KRT23", "VCAN", "FGF4", "DEFA5", "TANGO6",
              "FGF9", "FGF20", "DPEP1", "SEZ6L", "WNT5A", "WNT10A", "NKX2-5", "SLC1A5",
              "SEMA3C", "CAMK2B", "EFHD1", "KREMEN2", "LMO2", "NOTUM", "WNT10B", "TMEM45B",
              "DKK4", "RGPD3", "ITGA1", "PCCA", "CTNNB1", "MUC1")  # Wet epithelium genes

KE_genes <- c("SCUBE3", "TAGLN", "KIF5C", "CALB1", "ISG15", "CD24", "SYTL2", "CALML3",
              "ODAM", "HES1", "ATP6V0D2", "KRT75", "CLDN10", "MGP", "IFI6", "TM4SF1",
              "MACC1", "KRT6B", "VIM", "HHIP", "TRIB1", "PROCR", "GBP1", "MX1",
              "KRT16", "GJB6", "SLC20A2", "KLF6", "KRT6C", "IFIT1", "MAF", "FLNB",
              "GJB2", "CLDN1", "MYC", "KRT7")  # Keratinizing epithelium genes

PE_genes <- c("CTNNAL1", "SERTAD1", "CRABP1", "MT1E", "ATF3", "DLX6-AS1", "NFKBIA", "MT1X",
              "CYR61", "IFRD1", "EGR1", "MT1F", "NR4A1", "JUNB", "SLC39A10", "TSPAN8",
              "FOS", "MT1G", "KLF10", "CD47", "SFRP1", "FJX1", "EGR3", "EPHA7",
              "IRF1", "WIF1", "COL17A1", "FOSB", "FRZB", "IER3", "JUN", "JUND", "HSPA1A",
              "EGR2", "DDIT3", "POSTN", "ARL4D", "FAM41C", "SLC35G1", "KCTD12", "EPB41L2",
              "GADD45B", "PTCH1", "KLHL42", "VSTM2B")  # Pre-ameloblast epithelium genes

#---------------------------------------------
# Calculate Module Scores for Each Cell
#---------------------------------------------

# Join assay layers, overwrite existing if present
epi_cells_acp1_1 <- JoinLayers(epi_cells_acp1, overwrite = TRUE)

# Add module score for cell cycle genes
epi_cells_acp1_1 <- AddModuleScore(epi_cells_acp1_1, features = list(CC_genes), name = "CC_Score")

# Add module score for wet epithelium genes
epi_cells_acp1_1 <- AddModuleScore(epi_cells_acp1_1, features = list(WE_genes), name = "WE_Score")

# Add module score for keratinizing epithelium genes
epi_cells_acp1_1 <- AddModuleScore(epi_cells_acp1_1, features = list(KE_genes), name = "KE_Score")

# Add module score for pre-ameloblast epithelium genes
epi_cells_acp1_1 <- AddModuleScore(epi_cells_acp1_1, features = list(PE_genes), name = "PE_Score")

# View metadata before and after adding module scores
View(epi_cells_acp1@meta.data)
View(epi_cells_acp1_1@meta.data)

# Visualize module scores on UMAP, arranged in 2 columns, with title
epi_featureplot_acp1 <- FeaturePlot(
  epi_cells_acp1_1,
  features = c("CC_Score1", "WE_Score1", "KE_Score1", "PE_Score1"),
  ncol = 2
) + plot_annotation(title = "Acp1 Module Scores on Epithelial UMAP")

# Display module score plots
epi_featureplot_acp1

#---------------------------------------------
# Assign Dominant Epithelial Signature per Cell
#---------------------------------------------

# Define vector of score column names
score_cols <- c("CC_Score1", "WE_Score1", "KE_Score1", "PE_Score1")

# Extract scores from metadata
scores <- epi_cells_acp1_1@meta.data[, score_cols]

# For each cell, find index of the highest scoring signature
best_idx_acp1 <- apply(scores, 1, which.max)

# Assign the highest scoring signature's name to a new metadata column
epi_cells_acp1_1$epi_type <- score_cols[best_idx_acp1]

# Define colors for dominant epithelial types
epi_colors <- c(
  CC_Score1 = "lightblue",
  WE_Score1 = "darkgreen",
  KE_Score1 = "orange",
  PE_Score1 = "purple"
)

# Plot UMAP with cells colored by dominant epithelial signature
epi_umap_acp1 <- DimPlot(
  epi_cells_acp1_1,
  group.by = "epi_type",
  reduction = "umap",
  label = FALSE
) +
  ggtitle("Dominant Epithelial Signature per Cell") +
  scale_color_manual(values = epi_colors)

# Display the plot
epi_umap_acp1

# View metadata with new epi_type column
View(epi_cells_acp1_1@meta.data)

#---------------------------------------------
# Annotate Subclusters within Epithelial Cells
#---------------------------------------------

# Extract cluster IDs for epithelial subclusters at resolution 0.1
cluster_ids_epi_acp1 <- as.character(epi_cells_acp1$RNA_snn_res.0.1)

# Map cluster IDs to epithelial subtypes
celltype_annotation_acp1 <- dplyr::case_when(
  cluster_ids_epi_acp1 %in% c("1", "3") ~ "PE",
  cluster_ids_epi_acp1 == "0" ~ "KE",
  cluster_ids_epi_acp1 == "2" ~ "WE",
  cluster_ids_epi_acp1 == "4" ~ "CC"
)

# Add subcluster annotations to metadata
epi_cells_acp1$celltype <- celltype_annotation_acp1

# View metadata with subcluster annotations
View(epi_cells_acp1@meta.data)

# Plot UMAP colored by epithelial subtypes
Final_acp1 <- DimPlot(
  epi_cells_acp1,
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

#---------------------------------------------
# Save Raw Counts and Metadata
#---------------------------------------------

# Extract raw counts matrix from epithelial subset
counts_acp1 <- GetAssayData(epi_cells_acp1, assay = "RNA", slot = "counts")

# Convert counts matrix to data frame
counts_df_acp1 <- as.data.frame(as.matrix(counts_acp1))

# View raw counts matrix and data frame
View(counts_acp1)
View(counts_df_acp1)

# Extract metadata and add cell IDs as a column
metadata_df_acp1 <- epi_cells_acp1@meta.data
metadata_df_acp1$cell_id <- rownames(metadata_df_acp1)

# Save raw counts and metadata to CSV files in output directory
write.csv(counts_df_acp1, file = file.path(OUT_DIR, "acp1_epi_cells_raw_counts.csv"))
write.csv(metadata_df_acp1, file = file.path(OUT_DIR, "acp1_epi_cells_metadata.csv"), row.names = FALSE)

# Combine and display clustering and module-based annotation plots side-by-side
Plot_epi_cells_acp1 + epi_umap_acp1
Final_acp1 + epi_umap_acp1