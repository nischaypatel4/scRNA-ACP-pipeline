#---------------------------------------------
#---------------------------------------------
# Marker Expression in Integrated ACP-Tooth Object
#---------------------------------------------
#---------------------------------------------

#-------------------------------------------------------------------------------
# Libraries and Loading
#-------------------------------------------------------------------------------
library(Seurat)          # For single-cell RNA-seq analysis
library(dplyr)           # Data manipulation verbs
library(ggplot2)         # Visualization and plotting
library(patchwork)       # Combine multiple ggplots easily
library(Matrix)          # Sparse matrix operations
library(scales)          # Scaling functions for visualization
library(DropletUtils)    # Tools for single-cell RNA-seq data processing
library(harmony)         # Batch correction and data integration
library(RColorBrewer)    # Color palettes for plotting
library(biomaRt)         # Interface to BioMart databases
library(curl)            # HTTP requests utilities
library(AnnotationDbi)   # Annotation database interface
library(org.Hs.eg.db)    # Human gene annotation database
library(RColorBrewer)    # (Repeated) Color palettes for plotting

# Load integrated ACP Seurat object
acp_int <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_Integrated.rds")

# Load Tooth Seurat object
tooth <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Tooth.rds")

# Plot UMAP of ACP integrated object colored by 'subcelltype' with manual colors
Annotation_acp_final <- DimPlot(
  acp_int,
  group.by = "subcelltype",       # Group cells by subcelltype metadata
  reduction = "umap.harmony",     # Use Harmony corrected UMAP embeddings
  label = TRUE,                   # Add labels for clusters
  cols = c(                       # Manual color mapping for subcelltypes
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
Annotation_acp_final  # Display the plot

# View metadata of the ACP integrated object
View(acp_int@meta.data)

# 1) Create a new Seurat object subset containing only epithelial cell types specified
acp_epi <- subset(
  acp_int,
  subset = celltype %in% c("Epithelial_KE", "Epithelial_PE", "Epithelial_CC", "Epithelial_WE")
)

# Plot feature expression of oral epithelial markers (KRT6A, KRT13) on UMAP
OE <- FeaturePlot(acp_int, features = c("KRT6A", "KRT13"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Plot feature expression of enamel knot marker genes on UMAP
EK <- FeaturePlot(acp_int, features = c("SHH", "LEF1", "FGF4", "FGFR2", "CADM1"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Dental epithelial markers - commented out
# DE <- FeaturePlot(acp_int, features = c("PITX2", "KRT5", "KRT17"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Cervical Loop markers - commented out
# CL <- FeaturePlot(acp_int, features = c("LGR6", "TRPM3", "DACH1"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Outer enamel epithelium marker IGFBP5 expression plot
OEE <- FeaturePlot(acp_int, features = "IGFBP5", reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Inner enamel epithelium markers UNC5C and RYR2 expression plot
IEE <- FeaturePlot(acp_int, features = c("UNC5C", "RYR2"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Stellate intermediate odontoblast (SIO) marker - commented out
# SIO <- FeaturePlot(acp_int, features = c("CDH12"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Stellate reticulum intermediate (SRI) markers - commented out
# SRI <- FeaturePlot(acp_int, features = c("PCDH7", "FLNA", "SCUBE3"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Stellate intermediate (SII) markers FBN2, AFF3, AUTS2 expression plot
SII <- FeaturePlot(acp_int, features = c("FBN2", "AFF3", "AUTS2"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Stellate reticulum outer (SRO) marker ENOX1 expression plot
SRO <- FeaturePlot(acp_int, features = c("ENOX1"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Papillary layer (PA) markers VWDE, SHH, CYP2S1, KIF5C, RERE, MSI2 expression plot
PA <- FeaturePlot(acp_int, features = c("VWDE", "SHH", "CYP2S1", "KIF5C", "RERE", "MSI2"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Ameloblast (AM) markers DSPP, AMELX, AMBN expression plot
AM <- FeaturePlot(acp_int, features = c("DSPP", "AMELX", "AMBN"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Secretory ameloblast (SAM) markers MMP20, PLOD2, ENAM expression plot
SAM <- FeaturePlot(acp_int, features = c("MMP20", "PLOD2", "ENAM"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))


# Plot enamel knot markers on epithelial subset acp_epi UMAP
EK <- FeaturePlot(acp_epi, features = c("SHH", "LEF1", "FGF4", "FGFR2", "CADM1"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Plot OEE marker IGFBP5 on epithelial subset acp_epi UMAP
OEE <- FeaturePlot(acp_epi, features = "IGFBP5", reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Plot IEE markers UNC5C, RYR2 on epithelial subset acp_epi UMAP
IEE <- FeaturePlot(acp_epi, features = c("UNC5C", "RYR2"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Plot SII markers FBN2, AFF3, AUTS2 on epithelial subset acp_epi UMAP
SII <- FeaturePlot(acp_epi, features = c("FBN2", "AFF3", "AUTS2"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))

# Plot SRO marker ENOX1 on epithelial subset acp_epi UMAP
SRO <- FeaturePlot(acp_epi, features = c("ENOX1"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))


# Display the generated feature plots for enamel knot and epithelial subsets
EK
IEE
OEE
SII
SRO

## GO: 
# VCAN/ROBO1/ROBO2/CADM1 : Cell Recognition
# APCDD1/ROBO1/ROBO2 : regulation of animal organ morphogenesis
# ROBO1/ROBO2/CADM1 : homophilic cell adhesion via plasma membrane adhesion molecules

# Plot feature expression for genes involved in GO terms related to cell recognition and morphogenesis
FeaturePlot(acp_epi, features = c("VCAN", "CADM1", "ROBO1", "ROBO2", "APCDD1"), reduction = "umap.harmony", cols = c("lightgrey", "darkred"))