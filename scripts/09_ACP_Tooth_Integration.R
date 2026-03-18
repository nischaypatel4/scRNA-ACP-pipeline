#---------------------------------------------
#---------------------------------------------
# ACP and Tooth Integration
#---------------------------------------------
#---------------------------------------------

########################################################
# Libraries
######################################################## 
library(Seurat)         # For single-cell RNA-seq data analysis
library(dplyr)          # Data manipulation
library(ggplot2)        # Plotting
library(patchwork)      # Combine ggplot2 plots
library(Matrix)         # Sparse matrix operations
library(scales)         # Scale functions for visualization
library(DropletUtils)   # Single-cell data utilities
library(harmony)        # Batch correction integration
library(RColorBrewer)   # Color palettes for plots
library(biomaRt)        # BioMart database access for annotations
library(curl)           # Curl utilities for downloading data
library(AnnotationDbi)  # Database interface for annotation
library(org.Hs.eg.db)   # Human gene annotation database
library(RColorBrewer)   # (Repeated) Color palettes
library(GO.db)          # Gene Ontology database
library(clusterProfiler) # Functional enrichment analysis
library(ReactomePA)     # Reactome pathway analysis
library(org.Hs.eg.db)   # (Repeated) Human annotation database
library(enrichplot)     # Enrichment result visualization

proc_dir <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed"
# Set directory for processed data

########################################################
# ACP_Preprocessing
######################################################## 
acp_raw <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_Merged_Raw.rds")
# Load raw ACP dataset
acp_int <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_Integrated.rds")
# Load integrated ACP dataset

# Find intersecting cells between raw and integrated datasets
common_cells_acp <- intersect(Cells(acp_raw), Cells(acp_int))

# Initialize celltype column as NA for all cells in raw data
acp_raw$celltype <- NA

# Assign subcelltype values from integrated data to raw data for common cells
acp_raw$celltype[common_cells_acp] <- acp_int$subcelltype[common_cells_acp]

# Add prefix “acp_” to all entries in celltype column
acp_raw@meta.data$celltype <- paste0("acp_", acp_raw@meta.data$celltype)

View(acp_raw@meta.data)  # View metadata for confirmation

length(rownames(acp_raw)) # Count number of genes (rows) in raw ACP dataset

########################################################
# Tooth Preprocessing
######################################################## 
Tooth <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Tooth.rds")
# Load Tooth dataset

# Rename cluster_label column to celltype for consistency
colnames(Tooth@meta.data)[
  colnames(Tooth@meta.data) == "cluster_label"
] <- "celltype"

length(rownames(acp_raw)) # Confirm number of genes in ACP raw
length(rownames(Tooth))   # Confirm number of genes in Tooth data

# Find common genes between ACP and Tooth datasets
common_genes <- intersect(rownames(acp_raw), rownames(Tooth))
length(common_genes)  # Number of genes to keep after subsetting
any(duplicated(common_genes))  # Check for duplicate genes

# Subset ACP and Tooth datasets to common genes only
acp_raw <- subset(acp_raw, features = common_genes)
Tooth   <- subset(Tooth, features = common_genes)

length(rownames(acp_raw)) # Confirm gene count after subsetting
length(rownames(Tooth))   # Confirm gene count after subsetting

# Check unique cell types in Tooth dataset
unique(Tooth$celltype)
table(Tooth$celltype)
sort(table(Tooth$celltype), decreasing = TRUE)

# Subset Tooth data to exclude mesenchyme and salivary epithelium
Tooth1 <- subset(Tooth, subset = !celltype %in% c("Dental_mesenchyme", "Jaw_Mesenchyme", "Salivary_epithelium"))
View(Tooth1@meta.data)  # View Tooth metadata after filtering

########################################################
# Merging
######################################################## 
# Merge ACP raw and filtered Tooth data; add cell IDs to distinguish
merged <- merge(acp_raw, y = Tooth1, add.cell.ids = c("acp", "Tooth1"))

# Define ACP sample identifiers
acp_ids <- c("acp1", "acp2", "acp3")

# Create a new metadata column ‘identity’ based on orig.ident to distinguish ACP samples and Tooth cells
merged$identity <- ifelse(
  merged$orig.ident %in% acp_ids,
  merged$orig.ident,
  "tooth"
)

View(merged@meta.data)  # Inspect metadata after merging

########################################################
# SCTransform: Harmony
######################################################## 
set.seed(1234)  # Set seed for reproducibility

merged <- SCTransform(merged, verbose = TRUE)  # Normalize and variance-stabilize data using SCTransform
# saveRDS(merged, file = file.path(proc_dir, "ACP_Tooth_merged.rds"))  # Save processed object (commented out)

merged <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_Tooth_merged.rds")
# Load pre-saved SCTransform processed object

merged <- RunPCA(merged, verbose = TRUE)       # Compute principal components on normalized data

# Run UMAP dimensionality reduction on PCA embeddings
merged <- RunUMAP(merged, dims = 1:30, reduction.name = "umap")

# Plot UMAP colored by sample identity (unintegrated)
p0 <- DimPlot(merged, reduction = "umap", group.by = "identity",label=TRUE) + ggtitle("Unintegrated")

# Plot UMAP colored by celltype (unintegrated)
p1 <- DimPlot(merged, reduction = "umap", group.by = "celltype", label = TRUE) + ggtitle("Raw UMAP: celltype")

p0  # Display first plot
p1  # Display second plot

# Run Harmony integration on PCA embeddings to correct batch effects by sample identity
merged <- RunHarmony(
  object    = merged,
  group.by.vars = "identity",
  dims.use  = 1:30
)

# Run UMAP on Harmony corrected embeddings
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

# Plot Harmony UMAP colored by sample identity (integrated)
p2 <- DimPlot(merged, reduction = "umap.harmony", group.by = "identity",label=TRUE) + ggtitle("Integrated")

# Plot Harmony UMAP colored by celltype (integrated)
p3 <- DimPlot(merged, reduction = "umap.harmony", group.by = "celltype",label=TRUE) + ggtitle("Integrated")

p2  # Display integrated identity plot
p3  # Display integrated celltype plot

unique(merged$celltype)  # List unique celltypes after integration

# 1. Define celltypes to keep for focused plotting
keep_ct <- c(
  "acp_Epi_KE",  "acp_Epi_PE",  "acp_Epi_WE",  "acp_Epi_CC",
  "acp_Immune",  "Osteoblasts", "Otic_epithlium", "Oral_epithelium",
  "Dental_epithelium", "Immune_cells"
)

# 2. Subset merged object by celltypes of interest
to_plot <- subset(merged, subset = celltype %in% keep_ct)

# 3. Choose color palette with distinct colors for plotting
library(RColorBrewer)
palette <- brewer.pal(n = 12, name = "Paired")[1:length(keep_ct)]

# 4. Plot UMAP with selected celltypes colored by group
p1 <- DimPlot(
  to_plot,
  reduction  = "umap.harmony",
  group.by   = "celltype",
  cols       = palette,
  label      = FALSE,
  label.size = 4,
  pt.size    = 0.7
) + 
  ggtitle("UMAP: Selected ACP & Tooth Epithelial/Immune/Osteoblast Populations") +
  theme_minimal()

########################################################
# CCA_SCTransform
######################################################## 

# — 1) Split ACP raw object into individual samples by orig.ident (“acp1”, “acp2”, “acp3”)
acp_list  <- SplitObject(acp_raw,  split.by = "orig.ident")

# — 2) Split Tooth object into individual samples by identity (perio1, perio2, etc)
Tooth_list <- SplitObject(Tooth, split.by = "identity")

# — 3) Combine all individual datasets into one list
obj_list <- c(acp_list, Tooth_list)

# — 4) Run SCTransform normalization on each dataset separately
for (i in seq_along(obj_list)) {
  obj_list[[i]] <- SCTransform(obj_list[[i]], verbose = FALSE)
}

# — 5) Select integration features common to all datasets (top 3000 variable features)
features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = 3000)

# — 6) Prepare datasets for SCT integration with selected features
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)

# — 7) Run PCA on each dataset using selected features
obj_list <- lapply(obj_list, RunPCA, features = features, verbose = FALSE)

# — 8) Find integration anchors using reciprocal PCA (rpca) reduction method and SCT normalization
anchors <- FindIntegrationAnchors(
  object.list        = obj_list,
  normalization.method = "SCT",
  anchor.features      = features,
  reduction            = "rpca",
  dims                 = 1:30
)

# — 9) Integrate datasets using SCT normalization and identified anchors
cca_int <- IntegrateData(
  anchorset            = anchors,
  normalization.method = "SCT",
  dims                 = 1:30
)

# — 10) Set default assay to "integrated" for downstream analysis
DefaultAssay(cca_int) <- "integrated"

# — 11) Scale data, run PCA and UMAP on integrated assay
cca_int <- ScaleData(cca_int, verbose = FALSE)
cca_int <- RunPCA(cca_int, verbose = FALSE)
cca_int <- RunUMAP(cca_int, dims = 1:30)

# 1. Define celltypes of interest for plotting CCA integrated data
keep_ct <- c(
  "acp_Epi_KE",  "acp_Epi_PE",  "acp_Epi_WE",  "acp_Epi_CC",
  "acp_Immune",  "Osteoblasts", "Otic_epithlium", "Oral_epithelium",
  "Dental_epithelium", "Immune_cells"
)

# 2. Subset integrated object by selected celltypes
to_plot <- subset(cca_int, subset = celltype %in% keep_ct)

# 3. Use color palette from RColorBrewer
library(RColorBrewer)
palette <- brewer.pal(n = 12, name = "Paired")[1:length(keep_ct)]

# 4. Plot UMAP on integrated data
p2 <- DimPlot(
  to_plot,
  reduction  = "umap",
  group.by   = "celltype",
  cols       = palette,
  label      = FALSE,
  label.size = 4,
  pt.size    = 0.7
) + 
  ggtitle("UMAP: Selected ACP & Tooth Epithelial/Immune/Osteoblast Populations") +
  theme_minimal()

# Plot side-by-side plots for Harmony and CCA integration
p1 | p2

# Save Harmony and CCA integrated objects for future use
saveRDS(merged, file = file.path(proc_dir, "ACP_Tooth_Integrated_Harmony_SCT.rds"))
saveRDS(cca_int, file = file.path(proc_dir, "ACP_Tooth_Integrated_RPCA_SCT.rds"))

########################################################
# Read Merged object and plot
######################################################## 
merged <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/ACP_Tooth_Integrated_Harmony_SCT.rds")
# Load Harmony integrated merged object

View(merged@meta.data)  # View metadata

# 1. Read annotation CSV file containing new cluster labels
annotations <- read.csv("/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE184749/Tooth1_Epi_cell_annotations.csv", stringsAsFactors = FALSE)

# 2. Create named vector mapping cell IDs to new cluster labels
cluster_map <- setNames(annotations$new_cluster, annotations$cell)

# 3. Map new_cluster labels onto merged metadata using cell_id matching
merged@meta.data$new_cluster <- cluster_map[ merged@meta.data$cell_id ]

# 4. Check assignments (optional)
# table( !is.na(merged@meta.data$new_cluster) )
# head( merged@meta.data[, c("cell_id","new_cluster")] )
# View(merged@meta.data)

# Grab new_cluster and celltype columns as vectors for processing
nc <- merged@meta.data$new_cluster
ct <- merged@meta.data$celltype

# Find which entries are NA or empty string in new_cluster
empty_idx <- is.na(nc) | nc == ""

# Replace empty new_cluster entries with original celltype values
nc[empty_idx] <- ct[empty_idx]

# Update metadata column with corrected new_cluster values
merged@meta.data$new_cluster <- nc

# Sanity check (optional)
# table(empty_idx)               # Count of replaced entries
# head( merged@meta.data[, c("celltype","new_cluster")] )

# View final metadata table
# View(merged@meta.data)

# Plot UMAP grouping by celltype with labels and title
grouping <- DimPlot(merged, group.by = "celltype", reduction = "umap.harmony", label = TRUE) +
  ggtitle("Grouped Clusters with Manual Colors")

########################################################
# Clustering for FinallMarkers of that Small WE genes
######################################################## 

# Find nearest neighbors on Harmony embedding for clustering
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)

# Find clusters at multiple resolutions (0.1, 0.3, 0.5)
merged <- FindClusters(merged, resolution = c(0.1,0.3,0.5))

View(merged@meta.data)  # View metadata with cluster info

# Plot clusters at resolution 0.1
Clusters <- DimPlot(to_plot, group.by = "SCT_snn_res.0.1", label=TRUE, reduction = "umap.harmony") +
  ggtitle("Grouped Clusters with Manual Colors")

merged1 <- merged  # Copy merged object to merged1
View(merged1@meta.data)

# Subset merged1 to remove ACP sample identity labels
merged1_filtered <- subset(
  merged1,
  subset = !(identity %in% c("acp1", "acp2", "acp3"))
)

# Plot clusters of filtered data
Clusters_1 <- DimPlot(merged1_filtered, group.by = "SCT_snn_res.0.1", label=TRUE, reduction = "umap.harmony") +
  ggtitle("Grouped Clusters with Manual Colors")

# Set active identities to clustering results at resolution 0.1
Idents(merged1_filtered) <- "SCT_snn_res.0.1"

## total 14 clusters detected, next step: find markers for cluster 7 vs rest

# Run differential expression to find markers upregulated in cluster 7 compared to all others
cluster7_markers <- FindMarkers(
  merged1_filtered,
  ident.1 = "7",       # Target cluster to compare
  ident.2 = NULL,      # Compare against all other clusters
  only.pos = TRUE,     # Only positive markers (upregulated)
  min.pct = 0.25,      # Minimum percent of cells expressing gene
  logfc.threshold = 0.25
)

View(cluster7_markers)  # View marker genes for cluster 7

# Filter markers for statistical significance (adjusted p-value < 0.05)
cluster7_sig <- cluster7_markers[cluster7_markers$p_val_adj < 0.05, ]

# Extract gene names of significant markers
significant_genes <- rownames(cluster7_sig)

# View first few significant gene names
head(significant_genes)

# Convert gene names to list (optional for enrichment)
significant_genes_list <- as.list(significant_genes)

#-------------------------------------------------------------------------------
# GeneEnrichment
#-------------------------------------------------------------------------------
# 8. GO enrichment on the significant genes
ego_common <- enrichGO(
  gene          = significant_genes_list,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",             # Biological Process ontology
  pAdjustMethod = "BH",             # Benjamini-Hochberg p-value adjustment
  pvalueCutoff  = 0.05,             # P-value cutoff for significance
  qvalueCutoff  = 0.1               # Q-value cutoff for significance
)

# View enrichment results table
head(as.data.frame(ego_common))

# Plot enrichment results as dotplot showing top 20 categories
dotplot(ego_common, showCategory = 20) +
  ggtitle("GO BP Enrichment of Common DEGs")

View(as.data.frame(ego_common))  # View full enrichment results

# Example enriched pathways with genes:
# VCAN/ROBO1/ROBO2/CADM1 : Cell Recognition
# APCDD1/ROBO1/ROBO2 : regulation of animal organ morphogenesis
# ROBO1/ROBO2/CADM1 : homophilic cell adhesion via plasma membrane adhesion molecules

########################################################
# Plotting Merged object
######################################################## 

# View metadata for cluster names and counts
View(merged@meta.data)

# Sort unique clusters and print
sorted_clusters <- sort(unique(merged@meta.data$new_cluster))
print(sorted_clusters)

# Count and sort clusters by frequency
table(merged@meta.data$new_cluster)
cluster_counts <- sort(table(merged@meta.data$new_cluster), decreasing = TRUE)
print(cluster_counts)
View(cluster_counts)

library(Polychrome)  # Load Polychrome package for color palettes (if used)

View(merged@meta.data)  # View metadata again

# 1. Subset merged object by chosen cluster names
keep_ct <- c(
  "acp_Epi_KE","acp_Epi_PE","acp_Epi_WE",
  "acp_Epi_CC",
  "Ameloblasts",
  "IEE","OEE", "EK",
  "AM","PA","SAM",
  "SRI","SRO"
)
to_plot <- subset(merged, subset = new_cluster %in% keep_ct)

# Copy new_cluster to grouped_cluster column for manual grouping
to_plot$grouped_cluster <- to_plot$new_cluster

# Group related clusters into broader categories
to_plot$grouped_cluster[to_plot$new_cluster %in% c("Ameloblasts","AM","PA","SAM")] <- "Ameloblasts"
to_plot$grouped_cluster[to_plot$new_cluster %in% c("SRI","SRO")]              <- "Reticulum"
to_plot$grouped_cluster[to_plot$new_cluster %in% c("IEE")]             <- "Inner Enamel"
to_plot$grouped_cluster[to_plot$new_cluster %in% c("OEE")]             <- "Outer Enamel"
to_plot$grouped_cluster[to_plot$new_cluster %in% c("EK")]   <- "Enamel Knot"

# Manual color mapping for grouped clusters
cols <- c(
  "acp_Epi_PE" = "#0000FF",    # Dark Blue
  "acp_Epi_KE" = "#ff7f0e",    # Orange
  "acp_Epi_WE" = "#06402B",    # Dark Green
  "acp_Epi_CC" = "#ADD8E6",
  "Ameloblasts" = "#d62728",    # Red
  "Reticulum"  = "#000000",    # Black
  "Inner Enamel"= "#FFFF00",    # Yellow
  "Outer Enamel"= "#008000",    # Dark Green? (Note: original comment says Yellow but hex is green)
  "Enamel Knot"         = "#17becf"     # Grey (actually cyan/turquoise)
)

# Plot grouped clusters on UMAP with manual colors
grouping <- DimPlot(to_plot, group.by = "grouped_cluster", cols = cols, reduction = "umap.harmony") +
  ggtitle("Grouped Clusters with Manual Colors")

# Plot side-by-side: Clusters and grouping plots
Clusters | grouping

# Make sure active identity class is set to identity column
Idents(merged) <- merged$identity

########################################################
# Plotting ACP only object
######################################################## 
# 1) Subset merged object to only ACP cells (acp1, acp2, acp3)
acp_only <- subset(merged, idents = c("acp1","acp2","acp3"))

# 1. Subset merged to selected ACP epithelial clusters
keep_ct1 <- c(
  "acp_Epi_KE","acp_Epi_PE","acp_Epi_WE","acp_Epi_CC"
)
to_plot1 <- subset(merged, subset = new_cluster %in% keep_ct1)

# Subset ACP only to same clusters
to_plot2 <- subset(acp_only, subset = new_cluster %in% keep_ct1)

# Manual color mapping for ACP epithelial clusters
cols1 <- c(
  "acp_Epi_PE" = "#0000FF",    # Dark Blue
  "acp_Epi_KE" = "#ff7f0e",    # Orange
  "acp_Epi_WE" = "#06402B",    # Dark Green
  "acp_Epi_CC" = "#ADD8E6"
)

# Plot merged ACP clusters colored by new_cluster
acp_only_plot1 <- DimPlot(to_plot1, group.by = "new_cluster", cols = cols1, reduction = "umap.harmony") +
  ggtitle("Grouped Clusters with Manual Colors")

# Plot ACP only clusters colored by new_cluster
acp_only_plot2 <- DimPlot(to_plot2, group.by = "new_cluster", cols = cols1, reduction = "umap.harmony") +
  ggtitle("Grouped Clusters with Manual Colors")

# Combine plots side-by-side
acp_only_plot1 | acp_only_plot2

########################################################
# Plotting Tooth only 
######################################################## 
# 2) Subset merged to tooth cells only
tooth_only <- subset(merged, idents = "tooth")

# 1. Define tooth cluster groups to keep
keep_ct2 <- c(
  "Ameloblasts",
  "IEE","OEE", "EK",
  "AM","PA","SAM",
  "SRI","SRO"
)

# Subset merged object to these clusters
to_plot3 <- subset(merged, subset = new_cluster %in% keep_ct)

# Subset tooth only to same clusters
to_plot4 <- subset(tooth_only, subset = new_cluster %in% keep_ct)

# Copy new_cluster to grouped_cluster for manual grouping - merged object
to_plot3$grouped_cluster <- to_plot3$new_cluster
to_plot3$grouped_cluster[to_plot3$new_cluster %in% c("Ameloblasts","AM","PA","SAM")] <- "Ameloblasts"
to_plot3$grouped_cluster[to_plot3$new_cluster %in% c("SRI","SRO")]              <- "Reticulum"
to_plot3$grouped_cluster[to_plot3$new_cluster %in% c("IEE")]             <- "Inner Enamel"
to_plot3$grouped_cluster[to_plot3$new_cluster %in% c("OEE")]             <- "Outer Enamel"
to_plot3$grouped_cluster[to_plot3$new_cluster %in% c("EK")]   <- "Enamel Knot"

# Copy new_cluster to grouped_cluster for manual grouping - tooth only object
to_plot4$grouped_cluster <- to_plot4$new_cluster
to_plot4$grouped_cluster[to_plot4$new_cluster %in% c("Ameloblasts","AM","PA","SAM")] <- "Ameloblasts"
to_plot4$grouped_cluster[to_plot4$new_cluster %in% c("SRI","SRO")]              <- "Reticulum"
to_plot4$grouped_cluster[to_plot4$new_cluster %in% c("IEE")]             <- "Inner Enamel"
to_plot4$grouped_cluster[to_plot4$new_cluster %in% c("OEE")]             <- "Outer Enamel"
to_plot4$grouped_cluster[to_plot4$new_cluster %in% c("EK")]   <- "Enamel Knot"

# Manual color mapping for tooth clusters
cols2 <- c(
  "Ameloblasts" = "#d62728",    # Red
  "Reticulum"  = "#000000",    # Black
  "Inner Enamel"= "#FFFF00",    # Yellow
  "Outer Enamel"= "#008000",    # Dark Green
  "Enamel Knot"         = "#17becf"     # Grey (cyan)
)

# Plot merged tooth clusters
tooth_only_plot3 <- DimPlot(to_plot3, group.by = "grouped_cluster", cols = cols2, reduction = "umap.harmony") +
  ggtitle("Grouped Clusters with Manual Colors")

# Plot tooth-only clusters
tooth_only_plot4 <- DimPlot(to_plot4, group.by = "grouped_cluster", cols = cols2, reduction = "umap.harmony") +
  ggtitle("Grouped Clusters with Manual Colors")

# Display ACP plots side-by-side
acp_only_plot1 | acp_only_plot2

# Display tooth plots side-by-side
tooth_only_plot3 | tooth_only_plot4

# Display ACP plot next to tooth plot side-by-side
acp_only_plot1 | tooth_only_plot3 

# Stack ACP and tooth plots vertically
vertical_stack <- acp_only_plot1 / tooth_only_plot4
vertical_stack  # Display stacked plot

# Display tooth-only plots separately
tooth_only_plot3
tooth_only_plot4