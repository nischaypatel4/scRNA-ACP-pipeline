#---------------------------------------------
#---------------------------------------------
# ACP and Periondium Integration Visualisation
#---------------------------------------------
#---------------------------------------------

#-------------------------------------------------------------------------------
# Libraries and Loading
#-------------------------------------------------------------------------------

library(Seurat)          # Load Seurat package for single-cell analysis
library(dplyr)           # Load dplyr for data manipulation
library(ggplot2)         # Load ggplot2 for plotting
library(patchwork)       # Load patchwork to combine multiple ggplots
library(Matrix)          # Load Matrix package for sparse matrix handling
library(scales)          # Load scales for scaling functions in plots
library(DropletUtils)    # Load DropletUtils for single-cell data processing
library(harmony)         # Load harmony for batch correction and integration
library(RColorBrewer)    # Load RColorBrewer for color palettes

proc_dir <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed"  # Directory for processed data

# Load Seurat objects from disk
merged <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_PERIO_Integrated_Harmony1.rds")
cca_int_SCT <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_PERIO_Integrated_CCA_SCT1.rds")

# perio_int <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/PERIO_Integrated_1.rds")

# # 1) Prefix all subcelltype labels in perio_int
# perio_int$subcelltype <- paste0("perio_", perio_int$subcelltype)
# 
# # 2) Prefix all cell barcodes in perio_int
# #    This renames both the column names and the rownames of the meta.data
# new_barcodes <- paste0("perio_", colnames(perio_int))
# names(new_barcodes) <- colnames(perio_int)
# perio_int <- RenameCells(perio_int, new.names = new_barcodes)
# 
# # 3) Find cells common to both merged and perio_int
# common_cells <- intersect(colnames(merged), colnames(perio_int))
# 
# #    (Optional) Verify
# length(common_cells)
# # If you want to subset to only those cells:
# merged_common   <- subset(merged,   cells = common_cells)
# perio_common    <- subset(perio_int, cells = common_cells)
# 
# # 4) Transfer the (now‑prefixed) subcelltype into merged$celltype
# # Pull out the subcelltype vector from perio_int
# subtypes <- perio_int@meta.data[common_cells, "subcelltype", drop = TRUE]
# 
# # Overwrite merged$celltype for those shared cells
# merged@meta.data[common_cells, "celltype"] <- subtypes
# 
# # (Optional) Check
# table( merged@meta.data[common_cells, "celltype"] )
# 
# ####
# View(cca_int_SCT@meta.data)
# perio_int <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/PERIO_Integrated_1.rds")
# 
# # 1) Prefix all subcelltype labels in perio_int
# perio_int$subcelltype <- paste0("perio_", perio_int$subcelltype)
# # 3) Find cells common to both merged and perio_int
# common_cells <- intersect(colnames(cca_int_SCT), colnames(perio_int))
# 
# #    (Optional) Verify
# length(common_cells)
# # If you want to subset to only those cells:
# cca_int_SCT_common   <- subset(cca_int_SCT,   cells = common_cells)
# perio_common    <- subset(perio_int, cells = common_cells)
# 
# # 4) Transfer the (now‑prefixed) subcelltype into merged$celltype
# # Pull out the subcelltype vector from perio_int
# subtypes <- perio_int@meta.data[common_cells, "subcelltype", drop = TRUE]
# 
# # Overwrite merged$celltype for those shared cells
# cca_int_SCT@meta.data[common_cells, "celltype"] <- subtypes
# 
# # (Optional) Check
# table( cca_int_SCT@meta.data[common_cells, "celltype"] )
# ####

# saveRDS(merged, file =file.path(proc_dir, "ACP_PERIO_Integrated_Harmony1.rds"))
# saveRDS(cca_int_SCT, file = file.path(proc_dir, "ACP_PERIO_Integrated_CCA_SCT1.rds"))

set.seed(1234)            # Set random seed for reproducibility
View(perio_raw@meta.data) # View meta.data of perio_raw object (for inspection)
View(merged@meta.data)    # View meta.data of merged object (for inspection)

#-------------------------------------------------------------------------------
# Visualisation
#-------------------------------------------------------------------------------

#— 1) Define a robust palette for all categories in celltype & subcelltype
# Gather all unique levels of subcelltypes and celltypes from Harmony and CCA objects
h_sub    <- unique(merged$subcelltype)    # Harmony subcelltype (12)
h_cell   <- unique(merged$celltype)       # Harmony celltype (19)
c_sub    <- unique(cca_int_SCT$subcelltype)   # CCA subcelltype (13)
c_cell   <- unique(cca_int_SCT$celltype)      # CCA celltype (19)

# Pick a palette large enough for the max number (19) distinct categories:
bigpal   <- brewer.pal(12, "Set3")               # 12 colors from Set3 palette
bigpal   <- c(bigpal, brewer.pal(8, "Dark2"))    # Add 8 colors from Dark2 palette to make total 20 colors
# Truncate to exactly 19, or reorder as needed:
bigpal   <- bigpal[1:19]

# Function to map levels (categories) → colors
make_color_map <- function(levels_vec){
  # if fewer levels than palette, just take the first N colors
  pal <- bigpal[1:length(levels_vec)]
  names(pal) <- levels_vec
  pal
}

# Generate color maps for all categories
cols_h_sub  <- make_color_map(h_sub)    # Colors for Harmony subcelltype
cols_h_cell <- make_color_map(h_cell)   # Colors for Harmony celltype
cols_c_sub  <- make_color_map(c_sub)    # Colors for CCA subcelltype
cols_c_cell <- make_color_map(c_cell)   # Colors for CCA celltype

#— 2) UMAP by orig.ident / orig.ident1, Harmony vs CCA side‐by‐side

p_h_id  <- DimPlot(merged, reduction="umap.harmony", group.by="orig.ident1",  pt.size=0.4) + 
  ggtitle("Harmony UMAP — by orig.ident") + theme(plot.title=element_text(hjust=0.5))  # Plot Harmony UMAP colored by dataset origin
p_c_id  <- DimPlot(cca_int_SCT, reduction="umap",         group.by="orig.ident1", pt.size=0.4) + 
  ggtitle("CCA UMAP — by orig.ident1") + theme(plot.title=element_text(hjust=0.5))    # Plot CCA UMAP colored by dataset origin

(p_h_id | p_c_id) + plot_annotation(title="1) Dataset Separation: Harmony vs CCA")           # Combine plots side by side with annotation

#— 3) UMAP by subcelltype, Harmony vs CCA

p_h_sub  <- DimPlot(merged, reduction="umap.harmony", group.by="subcelltype",  cols=cols_h_sub,  pt.size=0.4) +
  ggtitle("Harmony — subcelltype") + theme(plot.title=element_text(hjust=0.5))             # Harmony UMAP colored by subcelltype
p_c_sub  <- DimPlot(cca_int_SCT, reduction="umap",         group.by="subcelltype",  cols=cols_c_sub,  pt.size=0.4) +
  ggtitle("CCA — subcelltype") + theme(plot.title=element_text(hjust=0.5))                 # CCA UMAP colored by subcelltype

(p_h_sub | p_c_sub) + plot_annotation(title="2) Subcelltype: Harmony vs CCA")               # Combine side by side with annotation

#— 4) Only Epithelial cells, then celltype in Harmony vs CCA

# define epithelial celltype labels in each object
# epi_h_cells <- c(
#   "acp_Epi_KE","acp_Epi_PE","acp_Epi_WE","acp_Epi_CC",
#   "perio_Epi1","perio_Epi2","perio_Epi3",
#   "perio_Epi4","perio_Epi5"
# )
# epi_c_cells <- c(
#   "acp_Epi_KE","acp_Epi_PE","acp_Epi_WE","acp_Epi_CC",
#   "perio_Epi1","perio_Epi2","perio_Epi3","perio_Epi4","perio_Epi5"
# )

# Define epithelial celltype labels used for subsetting
epi_h_cells <- c(
  "acp_Epi_KE","acp_Epi_PE","acp_Epi_WE","acp_Epi_CC",
  "perio_Epi1","perio_Epi2","perio_Epi3"
)
epi_c_cells <- c(
  "acp_Epi_KE","acp_Epi_PE","acp_Epi_WE","acp_Epi_CC",
  "perio_Epi1","perio_Epi2","perio_Epi3"
)

# subset by metadata column 'celltype' to get epithelial subsets
h_epi <- subset(merged,     subset = celltype %in% epi_h_cells)
c_epi <- subset(cca_int_SCT,    subset = celltype %in% epi_c_cells)

# UMAP by celltype in each subset
p_h_epi <- DimPlot(
  h_epi,
  reduction = "umap.harmony",
  group.by  = "celltype",
  cols      = cols_h_cell,
  pt.size   = 0.4
) + ggtitle("Harmony — Epithelial Only") +
  theme(plot.title = element_text(hjust=0.5))

p_c_epi <- DimPlot(
  c_epi,
  reduction = "umap",
  group.by  = "celltype",
  cols      = cols_c_cell,
  pt.size   = 0.4
) + ggtitle("CCA — Epithelial Only") +
  theme(plot.title = element_text(hjust=0.5))

# Combine epithelial plots side by side
(p_h_epi | p_c_epi) +
  plot_annotation(title = "3) Epithelial Subtypes: Harmony vs CCA")

# Feature plots of marker genes on Harmony UMAP
pa <- FeaturePlot(
  merged,
  features = "ODAM",
  reduction = "umap.harmony",
  order = TRUE
)
pa
pb <- FeaturePlot(
  merged,
  features = "WNT10A",
  reduction = "umap.harmony",
  order = TRUE
)
pb | p_h_epi

# 1. Define your gene list for cell cycle genes
CC_genes <- c(
  "KIAA0101","TOP2A","CDK1","UBE2C","NUSAP1","PBK","CENPM","TPX2",
  "BIRC5","ZWINT","CCNA2","ESCO2","CENPF","NCAPH","MKI67","CDCA3",
  "TK1","NDC80","CKAP2L","KIFC1","NCAPG","SGOL1","HMMR","NUF2",
  "AURKB","CEP55","CENPK","SPC25","CDKN3","CDCA5"
)

# 2. Add a module‐score for that gene set (here: “CC_Score”)
#    This computes, for each cell, a score = average expression of your genes
object <- AddModuleScore(
  object  = merged, 
  features = list(CC_genes), 
  name     = "CC_Score",
  assay    = "SCT",        # or whatever assay you normalized
  search   = FALSE         # set FALSE if your genes are already symbols
)

# The new meta.data column is called "CC_Score1"
head(object@meta.data[, "CC_Score1", drop=FALSE])

# 3. FeaturePlot of that combined score on Harmony UMAP
p <- FeaturePlot(
  object, 
  features       = "CC_Score1", 
  reduction      = "umap.harmony", 
  cols           = c("lightgrey", "darkred"), 
  pt.size        = 0.5
) + 
  scale_color_gradient(low = "lightgrey", high = "darkred") +
  ggtitle("Combined Cell‑Cycle Module Score")

print(p) | p_h_epi

#-------------------------------------------------------------------------------
# Cell Selection and Downstream Analysis: Cluster1
#-------------------------------------------------------------------------------

# There should a plot here
# UMAP by celltype in each subset
selection_plot <- DimPlot(
  h_epi,
  reduction = "umap.harmony",
  group.by  = "celltype",
  cols      = cols_h_cell,
  pt.size   = 0.4
) + ggtitle("Harmony — Epithelial Only") +
  theme(plot.title = element_text(hjust=0.5))

selection_plot

# Interactive cell selection: return vector of selected cells
selected_cells_1 <- CellSelector(plot = last_plot())

# Subset original object to only selected cells
selected_subset_1 <- subset(merged, cells = selected_cells_1)

# Add new metadata column Ident1 marking all cells as "selected"
selected_subset_1$Ident1 <- "selected"

View(selected_subset_1@meta.data)

# List of phellem markers to check
phellem_markers <- c("PER15", "PER49", "GPAT5", "CYP86A1", "HORST", 
                     "DAISY", "ASFT", "MYB67", "PBP1", "AT1G14120", "WOX4")

# Check if the genes are present in merged object
phellem_markers[phellem_markers %in% rownames(merged)]

# Feature plots of phellem markers on merged object
FeaturePlot(merged, features = phellem_markers, ncol = 3, order = TRUE)

FeaturePlot(merged, features = "WOX4", order = TRUE)

#-------------------------------------------------------------------------------
# Cell Selection and Downstream Analysis: Cluster2
#-------------------------------------------------------------------------------

# There should a plot here
# UMAP by celltype in each subset
selection_plot_1 <- DimPlot(
  h_epi,
  reduction = "umap.harmony",
  group.by  = "celltype",
  cols      = cols_h_cell,
  pt.size   = 0.4
) + ggtitle("Harmony — Epithelial Only") +
  theme(plot.title = element_text(hjust=0.5))

selection_plot_1

# Interactive cell selection: return vector of selected cells
selected_cells_2 <- CellSelector(plot = last_plot())

# Subset original object to only selected cells
selected_subset_2 <- subset(merged, cells = selected_cells_2)

# Add new metadata column Ident1 marking all cells as "selected"
selected_subset_2$Ident1 <- "selected"

View(selected_subset_1@meta.data)
View(selected_subset_2@meta.data)

#-------------------------------------------------------------------------------
# Processing
#-------------------------------------------------------------------------------

# Paths to input files for raw and integrated objects
acp_raw_path <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_Merged_Raw.rds"
acp_int_path <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_Integrated.rds"

perio_raw_path <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/PERIO_Merged_Raw.rds"
perio_int_path <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/PERIO_Integrated_1.rds"

# Load Seurat objects from the paths
acp_raw <- readRDS(acp_raw_path)
acp_int <- readRDS(acp_int_path)

perio_raw <- readRDS(perio_raw_path)
perio_int <- readRDS(perio_int_path)
#
# #-------------------------------------------------------------------------------
# # Preprocessing
# #-------------------------------------------------------------------------------
#
# Find intersecting cells between raw and integrated ACP objects
common_cells_acp <- intersect(Cells(acp_raw), Cells(acp_int))

# Initialize subcelltype column as NA for all cells in ACP raw object
acp_raw$celltype <- NA

# Add subcelltype values only for matching cells
acp_raw$celltype[common_cells_acp] <- acp_int$subcelltype[common_cells_acp]

# Find intersecting cells between raw and integrated PERIO objects
common_cells_perio <- intersect(Cells(perio_raw), Cells(perio_int))

# Initialize subcelltype column as NA for all cells in PERIO raw object
perio_raw$celltype <- NA

# Add subcelltype values only for matching cells
perio_raw$celltype[common_cells_perio] <- perio_int$subcelltype[common_cells_perio]

# Prefix subcelltype values with "acp_" or "perio_"
acp_raw@meta.data$celltype <- paste0("acp_", acp_raw@meta.data$celltype)
perio_raw@meta.data$celltype <- paste0("perio_", perio_raw@meta.data$celltype)

# View meta.data of raw ACP and PERIO objects
View(acp_raw@meta.data)
View(perio_raw@meta.data)

# Prefix ACP cell barcodes with "acp_"
acp_raw <- RenameCells(acp_raw, new.names = paste0("acp_", colnames(acp_raw)))

# Prefix PERIO cell barcodes with "perio_"
perio_raw <- RenameCells(perio_raw, new.names = paste0("perio_", colnames(perio_raw)))

# 2.1. Subset the raw Seurat objects to your selected cells based on cluster selections
cells_acp_Epi_KE   <- selected_subset_1$celltype == "acp_Epi_KE"
cells_perio_Epi1 <- selected_subset_1$celltype == "perio_Epi1"
cells_acp_Epi_PE   <- selected_subset_2$celltype == "acp_Epi_PE"
cells_perio_Epi2 <- selected_subset_2$celltype == "perio_Epi2"

# Get barcodes of the selected cells from each group
barcodes_acp_Epi_KE   <- Cells(selected_subset_1)[cells_acp_Epi_KE]
barcodes_perio_Epi1   <- Cells(selected_subset_1)[cells_perio_Epi1]
barcodes_acp_Epi_PE   <- Cells(selected_subset_2)[cells_acp_Epi_PE]
barcodes_perio_Epi2   <- Cells(selected_subset_2)[cells_perio_Epi2]

# Subset the raw Seurat objects based on selected barcodes
acp_Epi_KE_raw   <- subset(acp_raw, cells = barcodes_acp_Epi_KE)
acp_Epi_PE_raw   <- subset(acp_raw, cells = barcodes_acp_Epi_PE)
perio_Epi1_raw   <- subset(perio_raw, cells = barcodes_perio_Epi1)
perio_Epi2_raw   <- subset(perio_raw, cells = barcodes_perio_Epi2)

# View meta.data of each subsetted raw object
View(acp_Epi_KE_raw@meta.data)
View(acp_Epi_PE_raw@meta.data)
View(perio_Epi1_raw@meta.data)
View(perio_Epi2_raw@meta.data)

# Define marker gene lists for different oral epithelial subtypes
oral_markers <- c("KRT1", "KRT10", "IVL", "FLG", "LOR")
sulcular_markers <- c("KRT4", "KRT13", "DSG3", "CLDN1", "MUC1", "TGM3")
junctional_markers <- c("ODAM", "AMTN", "CXCL14", "ICAM1", "IL1B", "IL8", "CEACAM1", "MMP13")

# FeaturePlot for junctional markers on Harmony UMAP with red color gradient
FeaturePlot(merged, features = junctional_markers, reduction = "umap.harmony", 
            cols = c("lightgrey", "red"), ncol = 3) + 
  plot_annotation(title = "Oral Epithelium Markers")

# Save the raw subsetted objects to downstream folder for further analysis
saveRDS(acp_Epi_KE_raw, file =file.path("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Downstream", "acp_Epi_KE_raw.rds"))
saveRDS(acp_Epi_PE_raw, file =file.path("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Downstream", "acp_Epi_PE_raw.rds"))
saveRDS(perio_Epi1_raw, file =file.path("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Downstream", "perio_Epi1_raw.rds"))
saveRDS(perio_Epi2_raw, file =file.path("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Downstream", "perio_Epi2_raw.rds"))