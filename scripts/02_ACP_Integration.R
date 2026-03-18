#---------------------------------------------
#---------------------------------------------
# Integration of ACP 1,2,3
# Input: Final filtered Seurat objects (acp1_final, acp2_final, acp3_final)
# Output: ACP_Merged_Raw, ACP_Merged_SCT, ACP_Integrated
#---------------------------------------------
#---------------------------------------------

#---------------------------------------------
# Set a global seed for reproducibility (important for steps like RunUMAP, RunPCA)
set.seed(9867)
#---------------------------------------------

#---------------------------------------------
# Load required libraries for single-cell analysis and plotting
#---------------------------------------------
library(Seurat)         # Core single-cell toolkit
library(ggplot2)        # Data visualization
library(patchwork)      # Combining plots
library(sctransform)    # SCTransform normalization
library(harmony)        # Batch correction/integration
library(cowplot)        # Plot grid layouts
library(SeuratData)     # For SeuratData dataset management
library(patchwork)      # Loaded again (no harm, duplicates allowed)

#---------------------------------------------
# Set up file directories for processed data and outputs
#---------------------------------------------
PROC_DIR <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed"   # Path to processed data
OUT_DIR  <- "/conglilab/shared/projects/personal_projects/nischay25/project/results/ACP_Integration"   # Path for integration results
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)   # Create output directory if it doesn’t exist

#---------------------------------------------
# Load final filtered Seurat objects for ACP1, ACP2, ACP3
#---------------------------------------------
acp1 <- readRDS(file.path(PROC_DIR, "acp1_final.rds"))   # Load ACP1 object
acp2 <- readRDS(file.path(PROC_DIR, "acp2_final.rds"))   # Load ACP2 object
acp3 <- readRDS(file.path(PROC_DIR, "acp3_final.rds"))   # Load ACP3 object

#---------------------------------------------
# Merge the three objects into one Seurat object
#---------------------------------------------
combined <- merge(
  x = acp1,
  y = list(acp2, acp3),
  add.cell.ids = c("acp1", "acp2", "acp3"),   # Add dataset labels to cell names
  project = "acp_combined"                    # Name of the merged project
)

rm(acp1,acp2,acp3)   # Remove individual objects from memory to save RAM

saveRDS(combined, file = file.path(PROC_DIR, "ACP_Merged_Raw.rds"))   # Save merged raw object

#---------------------------------------------
# Standard Seurat preprocessing on merged object (before integration)
#---------------------------------------------
combined_swf <- NormalizeData(combined)                                         # Normalize counts
combined_swf <- FindVariableFeatures(combined_swf)                              # Find variable genes
combined_swf <- ScaleData(combined_swf)                                         # Scale data
combined_swf <- RunPCA(combined_swf, ndims = 50, verbose = TRUE, reduction.name = "pca")   # PCA
# ElbowPlot(combined_swf, ndims = 50,  reduction = "pca")                      # Optional: elbow plot (commented)
combined_swf <- RunUMAP(combined_swf, dims = 1:25, reduction = "pca", reduction.name = "umap.pca")   # UMAP
combined_swf <- FindNeighbors(combined_swf, dims = 1:25, reduction = "pca")     # Build neighbor graph
combined_swf <- FindClusters(combined_swf, resolution = c(0.1,0.2,0.3,0.5,1))   # Cluster cells at multiple resolutions

# Plot UMAP colored by sample origin (before integration)
p_unintegrated <- DimPlot(combined_swf, reduction = "umap.pca", group.by = "orig.ident") + 
  ggtitle("Unintegrated ACP Dataset") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave(file.path(OUT_DIR, "acp_unintegrated.png"), p_unintegrated, width = 10, height = 6)

p_unintegrated   # Display plot

#---------------------------------------------
# Integration using Harmony for batch correction
#---------------------------------------------
combined_HAR <- NormalizeData(combined)                                            # Normalize again
combined_HAR <- FindVariableFeatures(combined_HAR)                                 # Find variable features
combined_HAR <- ScaleData (combined_HAR)                                           # Scale data
combined_HAR <- RunPCA(combined_HAR, ndims = 50, verbose = TRUE, reduction.name = "pca")   # PCA again
# ElbowPlot(combined_HAR, ndims = 50,  reduction = "pca")                         # Optional elbow plot (commented)
combined_HAR <- RunHarmony(combined_HAR, group.by.vars = "orig.ident",reduction.name = "harmony")   # Run Harmony integration
harmony_embeddings <- Embeddings(combined_HAR, "harmony")                          # Extract Harmony embeddings
harmony_sds <- apply(harmony_embeddings, 2, sd)                                    # Compute standard deviations per Harmony dim
# plot(harmony_sds, type = "b", ...)                                              # Optional plot (commented)
combined_HAR <- RunUMAP(combined_HAR, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")   # UMAP on Harmony dims
combined_HAR <- FindNeighbors(combined_HAR, reduction = "harmony", dims = 1:30)    # Neighbors on Harmony dims
combined_HAR <- FindClusters(combined_HAR, resolution = c(0.1,0.2,0.3,0.5,1))      # Clustering

# Plot UMAP showing sample origins after integration
p_integrated <- DimPlot(combined_HAR, reduction = "umap.harmony", group.by = "orig.ident") + 
  ggtitle("Integrated ACP Dataset") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave(file.path(OUT_DIR, "acp_integrated.png"), p_integrated, width = 10, height = 6)

# Plot UMAP showing clusters after integration
p_cluster <- DimPlot(combined_HAR, reduction = "umap.harmony", group.by = "RNA_snn_res.0.1", label=TRUE) + 
  ggtitle("Clustered ACP Dataset") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
ggsave(file.path(OUT_DIR, "acp_clustered.png"), p_cluster, width = 10, height = 6)

# Display plots in R session
p_integrated
p_integrated + p_unintegrated
p_cluster

#---------------------------------------------
# Marker plots for major cell populations
#---------------------------------------------
# Generate FeaturePlots for each marker group
p_epi <- FeaturePlot(combined_HAR, features = c("CDH1", "EPCAM"), reduction = "umap.harmony", combine = FALSE)    # Epithelial
p_myelo <- FeaturePlot(combined_HAR, features = c("LYZ", "CD68"), reduction = "umap.harmony", combine = FALSE)   # Myeloid
p_t <- FeaturePlot(combined_HAR, features = c("CD3E"), reduction = "umap.harmony", combine = FALSE)              # T cells
p_b <- FeaturePlot(combined_HAR, features = c("CD79A"), reduction = "umap.harmony", combine = FALSE)             # B cells
p_astro <- FeaturePlot(combined_HAR, features = c("AQP4", "GFAP"), reduction = "umap.harmony", combine = FALSE)  # Astrocytes
p_opc <- FeaturePlot(combined_HAR, features = c("OLIG2"), reduction = "umap.harmony", combine = FALSE)           # OPC
p_oligo <- FeaturePlot(combined_HAR, features = c("MOG"), reduction = "umap.harmony", combine = FALSE)           # Oligodendrocytes
p_endo <- FeaturePlot(combined_HAR, features = c("PECAM1"), reduction = "umap.harmony", combine = FALSE)         # Endothelial

# Combine all FeaturePlots into a single grid
p_markers <- c(p_epi, p_myelo, p_t, p_b, p_astro, p_opc, p_oligo, p_endo)
p_markers <- plot_grid(plotlist = p_markers, ncol = 4)
ggsave(file.path(OUT_DIR, "acp_markers.png"), p_markers, width = 10, height = 6)

# Display combinations of cluster + marker plots
p_cluster + p_epi 
p_cluster + p_myelo 
p_cluster + p_t
p_cluster + p_b 
p_cluster + p_astro 
p_cluster + p_opc 
p_cluster + p_oligo 
p_cluster + p_endo 

#---------------------------------------------
# Cell Type Annotations
#---------------------------------------------
# Assign cluster IDs to broad cell types
cluster_ids_acp <- as.character(combined_HAR$RNA_snn_res.0.1)

celltype_annotation_acp <- dplyr::case_when(
  cluster_ids_acp %in% c("0", "1", "3", "5") ~ "Epithelial",
  cluster_ids_acp == "2" ~ "Myeloid",
  cluster_ids_acp == "7" ~ "Immune",
  cluster_ids_acp %in% c("4", "6") ~ "Astrocyte+OPC+Endothelial",
  cluster_ids_acp == "8" ~ "Oligodendrocyte",
  TRUE ~ "Unknown"
)

# Add annotated cell types into metadata
combined_HAR$celltype <- celltype_annotation_acp

# Plot annotated UMAP
Annotation_acp <- DimPlot(
  combined_HAR,
  group.by = "celltype",
  reduction = "umap.harmony",
  label = FALSE,
  cols = c(
    "Epithelial" = "blue",
    "Myeloid" = "orange",
    "Immune" = "red",
    "Astrocyte+OPC+Endothelial" = "pink",
    "Oligodendrocyte" = "brown",
    "Unknown" = "grey"
  )
)
Annotation_acp 
ggsave(file.path(OUT_DIR, "acp_celltype_annotation.png"), Annotation_acp , width = 10, height = 6)

Annotation_acp | p_cluster

#---------------------------------------------
# Epithelial AddModuleScore (4 epithelial signatures)
#---------------------------------------------
# Define 4 epithelial cell gene sets (CC, WE, KE, PE)
CC_genes <- c("KIAA0101","TOP2A","CDK1","UBE2C","NUSAP1","PBK","CENPM","TPX2",
              "BIRC5","ZWINT","CCNA2","ESCO2","CENPF","NCAPH","MKI67","CDCA3",
              "TK1","NDC80","CKAP2L","KIFC1","NCAPG","SGOL1","HMMR","NUF2",
              "AURKB","CEP55","CENPK","SPC25","CDKN3","CDCA5")

WE_genes <- c("MUCL1","RBMS3","COL1A2","KRT23","VCAN","FGF4","DEFA5","TANGO6",
              "FGF9","FGF20","DPEP1","SEZ6L","WNT5A","WNT10A","NKX2-5","SLC1A5",
              "SEMA3C","CAMK2B","EFHD1","KREMEN2","LMO2","NOTUM","WNT10B","TMEM45B",
              "DKK4","RGPD3","ITGA1","PCCA","CTNNB1","MUC1")

KE_genes <- c("SCUBE3","TAGLN","KIF5C","CALB1","ISG15","CD24","SYTL2","CALML3",
              "ODAM","HES1","ATP6V0D2","KRT75","CLDN10","MGP","IFI6","TM4SF1",
              "MACC1","KRT6B","VIM","HHIP","TRIB1","PROCR","GBP1","MX1",
              "KRT16","GJB6","SLC20A2","KLF6","KRT6C","IFIT1","MAF","FLNB",
              "GJB2","CLDN1","MYC","KRT7")

PE_genes <- c("CTNNAL1","SERTAD1","CRABP1","MT1E","ATF3","DLX6-AS1","NFKBIA","MT1X",
              "CYR61","IFRD1","EGR1","MT1F","NR4A1","JUNB","SLC39A10","TSPAN8",
              "FOS","MT1G","KLF10","CD47","SFRP1","FJX1","EGR3","EPHA7",
              "IRF1","WIF1","COL17A1","FOSB","FRZB","IER3","JUN","JUND","HSPA1A",
              "EGR2","DDIT3","POSTN","ARL4D","FAM41C","SLC35G1","KCTD12","EPB41L2",
              "GADD45B","PTCH1","KLHL42","VSTM2B")

# Combine gene sets into a list
epi_signatures <- list(
  Epithelial_CC = CC_genes,
  Epithelial_WE = WE_genes,
  Epithelial_KE = KE_genes,
  Epithelial_PE = PE_genes
)

# Subset epithelial cells for module scoring
epithelial <- subset(combined_HAR, subset = celltype == "Epithelial")
epithelial <- JoinLayers(epithelial, overwrite = TRUE)   # Join data layers

DefaultAssay(epithelial) <- "RNA"
epithelial <- AddModuleScore(
  object   = epithelial,
  features = epi_signatures,     # List of gene sets
  name     = "EpiModule",
  assay    = "RNA",
  layer    = "data",             # Joined data layer
  ctrl     = 100,
  seed.use = 42
)
head(epithelial@meta.data[, paste0("EpiModule", 1:4)])

# Rename EpiModule1..4 to CC, WE, KE, PE for clarity
module_cols <- paste0("EpiModule", seq_along(epi_signatures))
new_names   <- names(epi_signatures)
for (i in seq_along(module_cols)) {
  epithelial[[ new_names[i] ]] <- epithelial[[ module_cols[i] ]]
  epithelial[[ module_cols[i] ]] <- NULL
}

# Visualize module scores on UMAP
addmodule_plot <- FeaturePlot(
  epithelial,
  features = new_names,
  ncol     = 2
) +
  plot_annotation(title = "Epithelial module scores")
ggsave(file.path(OUT_DIR, "acp_epi_module_scores.png"), addmodule_plot , width = 10, height = 6)

#---------------------------------------------
# Epithelial Subcelltype Annotations
#---------------------------------------------
# Extract module score columns
score_cols <- new_names  # CC, WE, KE, PE
scores <- epithelial@meta.data[, score_cols]

# Find highest scoring signature per cell
best_idx <- apply(scores, 1, which.max)
epithelial$epi_type <- score_cols[best_idx]

# Assign colors for epithelial subtypes and plot UMAP
library(ggplot2)
epi_colors <- c(
  Epithelial_CC = "lightblue",
  Epithelial_WE = "darkgreen",
  Epithelial_KE = "orange",
  Epithelial_PE = "purple"
)

acp_epi_plot <- DimPlot(
  epithelial,
  group.by  = "epi_type",
  reduction = "umap.harmony",
  label     = TRUE
) + 
  ggtitle("Dominant Epithelial Signature per Cell") +
  scale_color_manual(values = epi_colors)
ggsave(file.path(OUT_DIR, "acp_epi_types.png"), acp_epi_plot , width = 10, height = 6)

# Identify epithelial cells in integrated object
epi_cells <- WhichCells(combined_HAR, expression = celltype == "Epithelial")

# Transfer epithelial subtypes into combined_HAR metadata
combined_HAR$celltype[epi_cells] <- epithelial$epi_type[epi_cells]

#---------------------------------------------
# Refine cluster-to-subcelltype mapping
#---------------------------------------------
cluster_ids_acp1 <- as.character(combined_HAR$RNA_snn_res.0.1)

celltype_annotation_acp1 <- dplyr::case_when(
  cluster_ids_acp1 == "0" ~ "Epi_PE",
  cluster_ids_acp1 == "1" ~ "Epi_KE",
  cluster_ids_acp1 == "3" ~ "Epi_WE",
  cluster_ids_acp1 == "5" ~ "Epi_CC",
  cluster_ids_acp1 == "2" ~ "Myeloid",
  cluster_ids_acp1 == "7" ~ "Immune",
  cluster_ids_acp1 %in% c("4", "6") ~ "Astrocyte+OPC+Endothelial",
  cluster_ids_acp1 == "8" ~ "Oligodendrocyte",
  TRUE ~ "Unknown"
)

# Add subcelltype labels
combined_HAR$subcelltype <- celltype_annotation_acp1

# combined_HAR <- readRDS("...")   # Commented out code for reloading object
# Plot final annotated UMAP
Annotation_acp_final <- DimPlot(
  combined_HAR,
  group.by = "subcelltype",
  reduction = "umap.harmony",
  cols = c(
    "Epi_PE" = "blue",
    "Epi_WE" = "darkgreen",
    "Epi_KE" = "orange",
    "Epi_CC" = "black",
    "Myeloid" = "yellow",
    "Immune" = "red",
    "Astrocyte+OPC+Endothelial" = "pink",
    "Oligodendrocyte" = "brown",
    "Unknown" = "grey"
  )
)
Annotation_acp_final
ggsave(file.path(OUT_DIR, "acp_subcelltype_annotation.png"), Annotation_acp_final , width = 10, height = 6)

#---------------------------------------------
# Save integrated Seurat object
#---------------------------------------------
saveRDS(combined_HAR, file = file.path(PROC_DIR, "ACP_Integrated.rds"))