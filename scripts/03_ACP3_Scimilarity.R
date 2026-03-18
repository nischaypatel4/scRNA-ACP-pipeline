#---------------------------------------------
#---------------------------------------------
# ACP Scimilarity Input Manipulation
# Input : acp3_final
# Output: acp3_epi_cells_raw_counts.csv and acp3_epi_cells_metadata.csv (Input to Scimilarity Tool)
#---------------------------------------------
#---------------------------------------------

#---------------------------------------------
# 0. Load libraries
#---------------------------------------------
library(Seurat)          # Load Seurat for single-cell analysis
library(ggplot2)         # Load ggplot2 for plotting
library(patchwork)       # Load patchwork for combining plots
library(sctransform)     # Load sctransform for normalization (though not used here explicitly)
library(harmony)         # Load harmony for batch correction (not used here explicitly)
library(cowplot)         # Load cowplot for plot grids

#---------------------------------------------
# Set directories
#---------------------------------------------
PROC_DIR <- "data/processed"  # Directory where processed data is stored
OUT_DIR  <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/Scimilarity_Input"  # Output directory for Scimilarity input
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)  # Create output directory if it doesn't exist

#---------------------------------------------
# 1. Load final filtered objects
#---------------------------------------------
acp3 <- readRDS(file.path(PROC_DIR, "acp3_final.rds"))  # Load the processed Seurat object

#---------------------------------------------
# acp3: Initial Preprocessing and Clustering
#---------------------------------------------

# Normalize raw counts
acp3_swf <- NormalizeData(acp3)  # Normalize the data (log-normalization by default)

# Identify variable features
acp3_swf <- FindVariableFeatures(acp3_swf)  # Find highly variable genes

# Scale data
acp3_swf <- ScaleData(acp3_swf)  # Scale and center features

# Perform PCA
acp3_swf <- RunPCA(acp3_swf, ndims = 50, verbose = FALSE)  # Run PCA for dimensionality reduction

# Plot Elbow Plot for PCA components
Elbow_acp3 <- ElbowPlot(acp3_swf, ndims = 50, reduction = "pca")  # Visualize variance explained by PCs
Elbow_acp3  # Display elbow plot

# Construct neighbor graph and cluster cells at multiple resolutions
acp3_swf <- FindNeighbors(acp3_swf, dims = 1:25, reduction = "pca")  # Find nearest neighbors using first 25 PCs
acp3_swf <- FindClusters(acp3_swf, resolution = c(0.1, 0.2, 0.3, 0.5, 0.7, 1))  # Cluster cells at multiple resolutions

# Run UMAP for visualization
acp3_swf <- RunUMAP(acp3_swf, dims = 1:25, reduction = "pca")  # Compute UMAP embeddings

# View metadata and set identities
View(acp3_swf@meta.data)  # Open metadata dataframe viewer
Idents(acp3_swf)          # Show current active identities
Idents(acp3_swf) <- "RNA_snn_res.0.2"  # Set cluster identities to resolution 0.2

# Plot UMAP colored by clusters
Plot_acp3 <- DimPlot(acp3_swf, group.by = c("RNA_snn_res.0.2"), ncol = 1, reduction = "umap")  # Plot clusters
Plot_acp3  # Display UMAP plot

#---------------------------------------------
# Marker Gene Expression Visualization
#---------------------------------------------

# Generate UMAP feature plots for known marker genes
p1_acp3 <- FeaturePlot(acp3_swf, features = c("CDH1", "EPCAM"), reduction = "umap", combine = FALSE)  # Epithelial markers
p2_acp3 <- FeaturePlot(acp3_swf, features = c("LYZ", "CD68"), reduction = "umap", combine = FALSE)   # Myeloid markers
p3_acp3 <- FeaturePlot(acp3_swf, features = c("CD3E"), reduction = "umap", combine = FALSE)          # T cell marker
p4_acp3 <- FeaturePlot(acp3_swf, features = c("CD79A"), reduction = "umap", combine = FALSE)         # B cell marker
p5_acp3 <- FeaturePlot(acp3_swf, features = c("AQP4", "GFAP"), reduction = "umap", combine = FALSE)  # Astrocyte markers
# p6_acp3 <- FeaturePlot(acp3_swf, features = c("SOX10", "OLIG2", "GFAP"), reduction = "umap", combine = FALSE)  # Commented out: Oligodendrocyte markers
p7_acp3 <- FeaturePlot(acp3_swf, features = c("PECAM1"), reduction = "umap", combine = FALSE)        # Endothelial marker

# Combine all marker plots into a grid
all_plots_acp3 <- c(p1_acp3, p2_acp3, p3_acp3, p4_acp3, p5_acp3, p7_acp3)  # Collect feature plots
combined_plot_acp3 <- plot_grid(plotlist = all_plots_acp3, ncol = 3)        # Arrange plots in grid with 3 columns

# Display combined plots and cluster UMAP
print(combined_plot_acp3)           # Show combined marker plots
combined_plot_acp3 | Plot_acp3     # Display side by side with cluster UMAP plot

#---------------------------------------------
# Annotate Cell Types from Clusters
#---------------------------------------------

# Map cluster IDs to known cell types
cluster_ids_acp3 <- as.character(acp3_swf$RNA_snn_res.0.2)  # Extract cluster IDs as characters

celltype_annotation_acp3 <- dplyr::case_when(    # Annotate clusters with cell type labels
  cluster_ids_acp3 %in% c("0", "1", "2", "3", "4", "6") ~ "Epithelial cells",  # Clusters assigned epithelial
  cluster_ids_acp3 == "5" ~ "Myeloid cells",                                   # Cluster 5 is myeloid
  cluster_ids_acp3 == "8" ~ "T cells",                                         # Cluster 8 is T cells
  cluster_ids_acp3 == "7" ~ "Astrocyte",                                       # Cluster 7 is astrocytes
  cluster_ids_acp3 == "9" ~ "Endothelial cells",                               # Cluster 9 endothelial
  TRUE ~ "Unknown"                                                             # All others unknown
)

# Add cell type labels to metadata
acp3_swf$celltype <- celltype_annotation_acp3  # Add annotation column to Seurat metadata

# Plot UMAP with annotated cell types
Annotation_acp3 <- DimPlot(
  acp3_swf,
  group.by = "celltype",     # Group by cell type annotation
  reduction = "umap",        # Use UMAP reduction
  label = FALSE,             # Do not label clusters
  cols = c(                  # Assign colors per cell type
    "Epithelial cells" = "blue",
    "Myeloid cells" = "orange",
    "T cells" = "red",
    "Astrocyte" = "pink",
    "Endothelial cells" = "black",
    "Unknown" = "grey"
  )
)
Annotation_acp3   # Show annotated UMAP plot
Annotation_acp3 | combined_plot_acp3  # Display side-by-side with marker plots

#---------------------------------------------
# Subset and Re-cluster Epithelial Cells
#---------------------------------------------

# Subset epithelial cells only
epi_cells_acp3 <- subset(acp3_swf, subset = celltype == "Epithelial cells")  # Keep only epithelial cells

# Perform downstream processing on epithelial subset
epi_cells_acp3 <- FindVariableFeatures(epi_cells_acp3)  # Find variable genes within epithelial subset
epi_cells_acp3 <- ScaleData(epi_cells_acp3)             # Scale the data
epi_cells_acp3 <- RunPCA(epi_cells_acp3, ndims = 50, verbose = FALSE)  # Run PCA on epithelial cells
ElbowPlot(epi_cells_acp3, ndims = 50, reduction = "pca")  # Elbow plot to determine PCs
epi_cells_acp3 <- FindNeighbors(epi_cells_acp3, dims = 1:25, reduction = "pca")  # Find neighbors on PCs 1-25
epi_cells_acp3 <- FindClusters(epi_cells_acp3, resolution = c(0.1, 0.2, 0.3, 0.5, 0.7, 1))  # Cluster epithelial cells
epi_cells_acp3 <- RunUMAP(epi_cells_acp3, dims = 1:25, reduction = "pca")  # Run UMAP for epithelial cells

# Set cluster identities at chosen resolution
Idents(epi_cells_acp3) <- "RNA_snn_res.0.1"  # Use clustering resolution 0.1 for epithelial cells
Plot_epi_cells_acp3 <- DimPlot(epi_cells_acp3, group.by = "RNA_snn_res.0.1", reduction = "umap", ncol = 1)  # Plot epithelial clusters
Plot_epi_cells_acp3 | Annotation_acp3  # Show epithelial clusters side-by-side with previous annotation

#---------------------------------------------
# Module Scores: Define Gene Signatures
#---------------------------------------------

# Define epithelial gene signatures
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

epi_cells_acp3_1 <- JoinLayers(epi_cells_acp3, overwrite = TRUE)  # Join layers, overwrite old if exists
epi_cells_acp3_1 <- AddModuleScore(epi_cells_acp3_1, features = list(CC_genes), name = "CC_Score")  # Add CC score
epi_cells_acp3_1 <- AddModuleScore(epi_cells_acp3_1, features = list(WE_genes), name = "WE_Score")  # Add WE score
epi_cells_acp3_1 <- AddModuleScore(epi_cells_acp3_1, features = list(KE_genes), name = "KE_Score")  # Add KE score
epi_cells_acp3_1 <- AddModuleScore(epi_cells_acp3_1, features = list(PE_genes), name = "PE_Score")  # Add PE score

View(epi_cells_acp3@meta.data)   # View original epithelial metadata
View(epi_cells_acp3_1@meta.data) # View metadata with module scores added

# Visualize module scores on UMAP
epi_featureplot_acp3 <- FeaturePlot(
  epi_cells_acp3_1,
  features = c("CC_Score1", "WE_Score1", "KE_Score1", "PE_Score1"),  # Use module score columns
  ncol = 2
) + plot_annotation(title = "acp3 Module Scores on Epithelial UMAP")  # Add title to combined plot
epi_featureplot_acp3  # Display module score plots

#---------------------------------------------
# Assign Dominant Epithelial Signature per Cell
#---------------------------------------------

# Identify the highest scoring signature per cell
score_cols <- c("CC_Score1", "WE_Score1", "KE_Score1", "PE_Score1")  # Score column names
scores <- epi_cells_acp3_1@meta.data[, score_cols]  # Extract score matrix
best_idx_acp3 <- apply(scores, 1, which.max)        # Find index of max score per cell
epi_cells_acp3_1$epi_type <- score_cols[best_idx_acp3]  # Assign highest scoring signature name

# Visualize dominant epithelial identity
epi_colors <- c(        # Define colors for signatures
  CC_Score1 = "lightblue",
  WE_Score1 = "darkgreen",
  KE_Score1 = "orange",
  PE_Score1 = "purple"
)

epi_umap_acp3 <- DimPlot(
  epi_cells_acp3_1,
  group.by = "epi_type",   # Group cells by dominant signature
  reduction = "umap",      # Use UMAP embedding
  label = FALSE            # No labels on plot
) +
  ggtitle("Dominant Epithelial Signature per Cell") +  # Add plot title
  scale_color_manual(values = epi_colors)               # Apply custom colors

epi_umap_acp3 | Plot_epi_cells_acp3  # Show side-by-side with epithelial clusters plot
View(epi_cells_acp3_1@meta.data)     # View metadata with epi_type annotation

#---------------------------------------------
# Annotate Subclusters within Epithelial Cells
#---------------------------------------------

# Map epithelial subcluster IDs to known epithelial types
cluster_ids_epi_acp3 <- as.character(epi_cells_acp3$RNA_snn_res.0.1)  # Get cluster IDs at res 0.1

celltype_annotation_acp3 <- dplyr::case_when(   # Assign cell type labels to clusters
  cluster_ids_epi_acp3 %in% c("0", "5") ~ "PE",  # Clusters 0 and 5 = PE
  cluster_ids_epi_acp3 %in% c("1", "2") ~ "KE",  # Clusters 1 and 2 = KE
  cluster_ids_epi_acp3 == "3" ~ "WE",             # Cluster 3 = WE
  cluster_ids_epi_acp3 == "4" ~ "CC"               # Cluster 4 = CC
)

# Add annotations
epi_cells_acp3$celltype <- celltype_annotation_acp3  # Add celltype info to metadata

View(epi_cells_acp3@meta.data)  # View updated metadata

# UMAP with epithelial annotations
Final_acp3 <- DimPlot(
  epi_cells_acp3,
  group.by = "celltype",     # Color by epithelial cell type annotation
  reduction = "umap",        # Use UMAP reduction
  label = FALSE,             # No labels
  cols = c(                  # Colors for each epithelial subtype
    "PE" = "blue",
    "KE" = "orange",
    "WE" = "darkgreen",
    "CC" = "pink",
    "Unknown" = "grey"
  )
)
Final_acp3  # Show final epithelial annotation UMAP plot

#---------------------------------------------
# Save Raw Counts and Metadata
#---------------------------------------------

# Extract raw counts matrix
counts_acp3 <- GetAssayData(epi_cells_acp3, assay = "RNA", slot = "counts")  # Get raw counts
counts_df_acp3 <- as.data.frame(as.matrix(counts_acp3))                       # Convert to data frame

View(counts_acp3)       # View raw counts matrix
View(counts_df_acp3)    # View counts dataframe

# Extract and update metadata
metadata_df_acp3 <- epi_cells_acp3@meta.data  # Extract metadata dataframe
metadata_df_acp3$cell_id <- rownames(metadata_df_acp3)  # Add cell IDs as a column

# Save to CSV
write.csv(counts_df_acp3, file = file.path(OUT_DIR, "acp3_epi_cells_raw_counts.csv"))  # Save counts CSV
write.csv(metadata_df_acp3, file = file.path(OUT_DIR, "acp3_epi_cells_metadata.csv"), row.names = FALSE)  # Save metadata CSV

# Combine clustering and module-based annotation plots
Plot_epi_cells_acp3 + epi_umap_acp3  # Plot epithelial clusters with module scores side-by-side
Final_acp3 + epi_umap_acp3            # Plot final annotation with module scores side-by-side