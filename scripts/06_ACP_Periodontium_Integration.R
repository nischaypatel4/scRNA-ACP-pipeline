#---------------------------------------------
#---------------------------------------------
# ACP and Periondium Integration using Harmony and Seurat's CCA
# Input : ACP and Periondium Merged Objects
# Output: Integrated Objects
#---------------------------------------------
#---------------------------------------------

#-------------------------------------------------------------------------------
# Libraries and Loading
#-------------------------------------------------------------------------------

library(Seurat)        # Single-cell RNA-seq analysis
library(dplyr)         # Data manipulation
library(ggplot2)       # Visualization
library(patchwork)     # Combine multiple ggplots
library(Matrix)        # Matrix operations
library(scales)        # Plot scales
library(DropletUtils)  # Processing droplet-based scRNA-seq data
library(harmony)       # Batch effect correction
library(RColorBrewer)  # Color palettes

proc_dir <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed"  # Directory for processed data

# Paths to input files
acp_raw_path <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_Merged_Raw.rds"   # ACP raw merged Seurat object
acp_int_path <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_Integrated.rds"   # ACP integrated Seurat object

perio_raw_path <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/PERIO_Merged_Raw.rds"  # Periodontium raw merged Seurat object
perio_int_path <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/PERIO_Integrated.rds"  # Periodontium integrated Seurat object

# Load Seurat objects
acp_raw <- readRDS(acp_raw_path)   # Load raw ACP object
acp_int <- readRDS(acp_int_path)   # Load integrated ACP object

perio_raw <- readRDS(perio_raw_path)  # Load raw Periodontium object
perio_int <- readRDS(perio_int_path)  # Load integrated Periodontium object

#-------------------------------------------------------------------------------
# Preprocessing
#-------------------------------------------------------------------------------

# Find intersecting cells between raw and integrated ACP objects
common_cells_acp <- intersect(Cells(acp_raw), Cells(acp_int))

# Initialize subcelltype column as NA for all ACP raw cells
acp_raw$celltype <- NA

# Add subcelltype values to raw ACP metadata only for matching cells
acp_raw$celltype[common_cells_acp] <- acp_int$subcelltype[common_cells_acp]

# Find intersecting cells between raw and integrated Periodontium objects
common_cells_perio <- intersect(Cells(perio_raw), Cells(perio_int))

# Initialize subcelltype column as NA for all Periodontium raw cells
perio_raw$celltype <- NA

# Add subcelltype values to raw Periodontium metadata only for matching cells
perio_raw$celltype[common_cells_perio] <- perio_int$subcelltype[common_cells_perio]

# Prefix celltype values with dataset label for ACP and Periodontium
acp_raw@meta.data$celltype <- paste0("acp_", acp_raw@meta.data$celltype)
perio_raw@meta.data$celltype <- paste0("perio_", perio_raw@meta.data$celltype)

# Merge both datasets
# merged <- merge(acp_raw, y = perio_raw, add.cell.ids = c("acp", "perio"))

# # Rename in 'merged'
# merged$celltype <- gsub(
#   pattern     = "^perio_Epithelial(\\d+)$",
#   replacement = "perio_Epi\\1",
#   x           = merged$celltype
# )

# merged$orig.ident1 <- ifelse(
#   merged$orig.ident %in% c("acp1", "acp2", "acp3"),
#   "acp",
#   ifelse(
#     merged$orig.ident %in% c("perio1", "perio2", "perio3", "perio4", "perio5"),
#     "perio",
#     merged$orig.ident  # fallback in case other identities exist
#   )
# )
# 
# merged$subcelltype <- ifelse(
#   merged$celltype %in% c("acp_Epi_KE", "acp_Epi_CC", "acp_Epi_WE", "acp_Epi_PE"),
#   "acp_EPI",
#   ifelse(
#     merged$celltype %in% c("perio_Epithelial1", "perio_Epithelial2", "perio_Epithelial3", "perio_Epithelial4", "perio_Epithelial5"),
#     "perio_EPI",
#     merged$celltype  # keep others unchanged
#   )
# )

#-------------------------------------------------------------------------------
# Pre‑integration UMAPs (for comparison)
#-------------------------------------------------------------------------------

# merged <- SCTransform(merged, verbose = FALSE)  # normalize (commented out here)
# saveRDS(merged, file = file.path(proc_dir, "ACP_PERIO_Merged_SCT.rds"))  # save merged object after normalization

set.seed(1234)  # Set seed for reproducibility
merged <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_PERIO_Merged_SCT.rds")  # Load SCT-normalized merged object

# Rename celltype to unify perio epithelial numbering
merged$celltype <- gsub(
  pattern     = "^perio_Epithelial(\\d+)$",
  replacement = "perio_Epi\\1",
  x           = merged$celltype
)

View(merged@meta.data)  # View metadata table

# # Check the result
# unique(merged$celltype)
# # unique(cca_int$celltype)

merged <- RunPCA(merged, verbose = TRUE)       # Perform PCA on merged data
# UMAP on raw PCs
merged <- RunUMAP(merged, dims = 1:30, reduction.name = "umap.raw")  # UMAP using first 30 PCs

# Plot UMAPs by different groupings
p0 <- DimPlot(merged, reduction = "umap.raw", group.by = "orig.ident1") + ggtitle("Unintegrated")
p1 <- DimPlot(merged, reduction = "umap.raw", group.by = "orig.ident") + ggtitle("Raw UMAP: Dataset")
p2 <- DimPlot(merged, reduction = "umap.raw", group.by = "celltype", label = FALSE) + ggtitle("Raw UMAP: celltype")
p3 <- DimPlot(merged, reduction = "umap.raw", group.by = "subcelltype", label = TRUE) + ggtitle("Raw UMAP: subcelltype")

#-------------------------------------------------------------------------------
# Harmony integration
#-------------------------------------------------------------------------------
# run Harmony on the SCTransformed PCA
merged <- RunHarmony(
  object    = merged,
  group.by.vars = "orig.ident",      # Batch variable to correct
  dims.use  = 1:30                   # Use first 30 PCs
)
# compute UMAP on the Harmony embeddings
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
# merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)   # neighbor graph (commented out)
# merged <- FindClusters(merged, resolution = 0.5)                     # clustering (commented out)

saveRDS(merged, file = file.path(proc_dir, "ACP_PERIO_Integrated_Harmony.rds"))  # Save Harmony integrated object

#-------------------------------------------------------------------------------
# CCA integration_LogNormalisation
#-------------------------------------------------------------------------------

# — 1) Split each merged object into individual samples
acp_list  <- SplitObject(acp_raw,  split.by = "orig.ident")   # Split ACP raw object by sample: acp1, acp2, acp3
perio_list <- SplitObject(perio_raw, split.by = "orig.ident")  # Split perio raw object by sample: perio1 to perio5

# — 2) Put them into one list
obj_list <- c(acp_list, perio_list)   # Combine both lists into one

# — 3) Standard log‐normalization & HVG selection on each
for (i in seq_along(obj_list)) {
  obj_list[[i]] <- NormalizeData(obj_list[[i]], normalization.method = "LogNormalize", scale.factor = 1e4)  # Normalize counts
  obj_list[[i]] <- FindVariableFeatures(obj_list[[i]], selection.method = "vst", nfeatures = 2000)        # Find HVGs
}

# — 4) Select integration features (genes) across all objects
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)  # Select 3000 integration features

# — 5) Run scaling & PCA on each object (needed for CCA)
for (i in seq_along(obj_list)) {
  obj_list[[i]] <- ScaleData(obj_list[[i]], features = features, verbose = FALSE)  # Scale selected features
  obj_list[[i]] <- RunPCA(obj_list[[i]], features = features, verbose = FALSE)    # PCA on selected features
}

# — 6) Find CCA‐based integration anchors
anchors <- FindIntegrationAnchors(
  object.list = obj_list,
  normalization.method = "LogNormalize",  # Match normalization used above
  anchor.features      = features,
  reduction            = "cca",           # Use canonical correlation analysis for anchors
  dims                 = 1:30
)

# Inspect a few anchors to see which cells and genes are being matched:
head(anchors@anchors)         # Cell-pair indices for anchors
head(anchors@anchor.features) # Top genes driving anchor identification

# — 7) Integrate data
cca_int <- IntegrateData(
  anchorset            = anchors,
  normalization.method = "LogNormalize",
  dims                 = 1:30
)

# — 8) Downstream on the integrated assay
DefaultAssay(cca_int) <- "integrated"  # Set integrated assay as default
cca_int <- ScaleData(cca_int, verbose = FALSE)  # Scale integrated data
cca_int <- RunPCA(cca_int, verbose = FALSE)     # PCA on integrated data
cca_int <- RunUMAP(cca_int, dims = 1:30)        # UMAP embedding

# Add simplified orig.ident1 grouping (ACP vs Periodontium)
cca_int$orig.ident1 <- ifelse(
  cca_int$orig.ident %in% c("acp1", "acp2", "acp3"),
  "acp",
  ifelse(
    cca_int$orig.ident %in% c("perio1", "perio2", "perio3", "perio4", "perio5"),
    "perio",
    cca_int$orig.ident  # fallback in case other identities exist
  )
)

# Simplify celltype to broader subcelltype groups
cca_int$subcelltype <- ifelse(
  cca_int$celltype %in% c("acp_Epi_KE", "acp_Epi_CC", "acp_Epi_WE", "acp_Epi_PE"),
  "acp_EPI",
  ifelse(
    cca_int$celltype %in% c("perio_Epi1", "perio_Epi2", "perio_Epi3", "perio_Epi4", "perio_Epi5"),
    "perio_EPI",
    cca_int$celltype  # keep others unchanged
  )
)

View(cca_int@meta.data)  # View metadata
saveRDS(anchors, file = file.path(proc_dir, "cca_anchors_lognorm.rds"))   # Save anchors
saveRDS(cca_int, file = file.path(proc_dir, "ACP_PERIO_Integrated_CCA.rds"))  # Save integrated CCA object

#-------------------------------------------------------------------------------
# CCA integration_SCTransform
#-------------------------------------------------------------------------------

# — 1) Split each merged object into individual samples
acp_list  <- SplitObject(acp_raw,  split.by = "orig.ident")   # Split ACP raw object by sample
perio_list <- SplitObject(perio_raw, split.by = "orig.ident") # Split Periodontium raw object by sample

# — 2) Put them into one list
obj_list <- c(acp_list, perio_list)  # Combine all objects

# — 3) Run SCTransform on each dataset (regression optional)
for (i in seq_along(obj_list)) {
  obj_list[[i]] <- SCTransform(obj_list[[i]], verbose = FALSE)  # SCTransform normalization
}

# — 4) Select integration features
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)  # Features for integration

# — 5) Prepare for SCT integration
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)  # Prepare data for SCT integration

# — 6) Run PCA (SCT stores it in the corrected assay)
obj_list <- lapply(obj_list, RunPCA, features = features, verbose = FALSE)  # PCA on SCT assay

# — 7) Find integration anchors using CCA
anchors <- FindIntegrationAnchors(
  object.list        = obj_list,
  normalization.method = "SCT",           # SCT normalization method
  anchor.features      = features,
  reduction            = "cca",           # Use CCA reduction for anchors
  dims                 = 1:30
)

# — 8) Integrate data
cca_int <- IntegrateData(
  anchorset            = anchors,
  normalization.method = "SCT",   # Must match above normalization
  dims                 = 1:30
)

# — 9) Downstream analysis on the "integrated" assay
DefaultAssay(cca_int) <- "integrated"   # Set integrated assay as default
cca_int <- ScaleData(cca_int, verbose = FALSE)    # Scale integrated data
cca_int <- RunPCA(cca_int, verbose = FALSE)       # PCA
cca_int <- RunUMAP(cca_int, dims = 1:30)          # UMAP

# Assign simplified orig.ident1 grouping
cca_int$orig.ident1 <- ifelse(
  cca_int$orig.ident %in% c("acp1", "acp2", "acp3"),
  "acp",
  ifelse(
    cca_int$orig.ident %in% c("perio1", "perio2", "perio3", "perio4", "perio5"),
    "perio",
    cca_int$orig.ident  # fallback for other identities
  )
)

# Assign simplified subcelltype groups
cca_int$subcelltype <- ifelse(
  cca_int$celltype %in% c("acp_Epi_KE", "acp_Epi_CC", "acp_Epi_WE", "acp_Epi_PE"),
  "acp_EPI",
  ifelse(
    cca_int$celltype %in% c("perio_Epi1", "perio_Epi2", "perio_Epi3", "perio_Epi4", "perio_Epi5"),
    "perio_EPI",
    cca_int$celltype  # Keep others unchanged
  )
)

View(cca_int@meta.data)   # View metadata after integration
saveRDS(anchors, file = file.path(proc_dir, "cca_anchors_SCT.rds"))          # Save SCT anchors
saveRDS(cca_int, file = file.path(proc_dir, "ACP_PERIO_Integrated_SCT.rds"))  # Save SCT integrated object