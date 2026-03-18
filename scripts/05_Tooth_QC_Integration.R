#---------------------------------------------
#---------------------------------------------
# Tooth Developmental Dataset QC,Filtering and Integration
# Input : Tooth Dataset
# Output: Processed Tooth Dataset
#---------------------------------------------
#---------------------------------------------

########################################################
# Libraries
######################################################## 
library(Seurat)        # Single-cell RNA-seq data analysis toolkit
library(Matrix)        # Sparse and dense matrix classes and methods
library(dplyr)         # Data manipulation and piping
library(ggplot2)       # Plotting system for R
library(patchwork)     # Combine ggplots easily
library(Matrix)        # Loaded again (redundant but harmless)
library(scales)        # Scale functions for visualization
library(DropletUtils)  # Utilities for droplet-based scRNA-seq data
library(harmony)       # Batch correction and integration
########################################################
proc_dir  <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed"
########################################################
# Loading Data
######################################################## 
# 1. Load everything
count_matrix <- Matrix::readMM("/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE184749/GSE184749_all_data_counts.mtx.gz")  # Read sparse matrix counts
gene_anno <- read.csv("/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE184749/GSE184749_all_data_gene_annotations.csv.gz")  # Gene annotations
cell_anno <- read.csv("/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE184749/GSE184749_all_data_cell_annotations.csv.gz")  # Cell annotations

rownames(gene_anno) <- gene_anno[,1]    # Assign rownames as gene IDs
row.names(cell_anno) <- cell_anno[,1]   # Assign rownames as cell IDs

length(rownames(gene_anno))  # Number of genes
length(rownames(cell_anno))  # Number of cells

# 2. Assign names to the matrix
rownames(count_matrix) <- rownames(gene_anno)    # Name matrix rows as genes
colnames(count_matrix) <- rownames(cell_anno)    # Name matrix columns as cells

# 3. (Optional) Convert to dgCMatrix explicitly to suppress coercion warning
count_matrix <- as(count_matrix, "dgCMatrix")

# 4. Re-create your Seurat object
Tooth <- CreateSeuratObject(
  counts    = count_matrix,      # Raw counts matrix
  project   = "Tooth",           # Project name
  meta.data = cell_anno          # Cell metadata
)

# 5. Quick re-check
head(rownames(Tooth))            # Check gene IDs
head(colnames(Tooth))            # Check cell barcodes
Assays(Tooth)                   # Check assays present ("RNA")
Layers(Tooth)                   # Check assay layers ("counts")
dim(GetAssayData(Tooth, assay="RNA", layer="counts"))  # Genes × cells dimensions
nrow(Tooth@meta.data) == ncol(GetAssayData(Tooth, "RNA","counts"))  # TRUE, confirm metadata matches counts
View(Tooth@meta.data)            # View metadata

# Filter out unwanted clusters and fix typo in cell type annotation
Tooth_subset <- {
  tmp <- subset(
    x      = Tooth,
    subset = !(cluster_number %in% c(5, 6, 7, 14, 18, 19, 20))  # Remove specific clusters
  )
  tmp$cluster_label[tmp$cluster_label == "Ednothelial_cells"] <- "Endothelial_cells"  # Correct typo in cluster label
  tmp$identity <- "tooth"        # Add new column with value “tooth” for all cells
  tmp
}

# Count how many rownames look like Ensembl IDs
ensembl_like <- grepl("^ENSG", rownames(Tooth_subset))    # Check if gene names start with ENSG
prop_ensembl <- sum(ensembl_like) / length(ensembl_like)  # Calculate proportion of Ensembl IDs

if (prop_ensembl > 0.9) {
  message("Tooth Seurat object is using Ensembl Gene IDs.")  # Mostly Ensembl IDs used
} else {
  message("Tooth Seurat object is using gene symbols.")      # Otherwise gene symbols
}

# 1. Get sparse counts and map IDs
counts_mat <- GetAssayData(Tooth_subset, assay = "RNA", slot = "counts")  # Extract raw counts matrix as dgCMatrix
ensembl_ids <- rownames(counts_mat)                                        # Ensembl IDs

# Map Ensembl IDs to gene symbols using org.Hs.eg.db annotation package
symbol_map <- mapIds(org.Hs.eg.db,
                     keys    = ensembl_ids,
                     column  = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

# 2. Drop unmapped and rename rows
keep    <- !is.na(symbol_map)              # Keep only genes with mapped symbols
cm_sub  <- counts_mat[keep, ]              # Subset count matrix to kept genes
groups  <- symbol_map[keep]                 # Corresponding gene symbols
rownames(cm_sub) <- groups                  # Rename rows to gene symbols

# 3. Collapse duplicates via rowsum (aggregate counts of duplicated gene symbols)
collapsed_mat <- rowsum(cm_sub, group = groups)

# Create new assay object with collapsed matrix and assign to Seurat object
new_rna_assay <- CreateAssayObject(counts = collapsed_mat)
Tooth_subset[["RNA"]] <- new_rna_assay

rownames(Tooth_subset)                     # Confirm updated rownames
Idents(Tooth_subset) <- Tooth_subset$identity  # Set identities to “tooth” for all cells

########################################################
# SCTransform
######################################################## 
# 2️⃣ SCTransform (replaces NormalizeData + FindVariableFeatures + ScaleData)
Tooth_subset <- SCTransform(
  object  = Tooth_subset,
  assay   = "RNA",
  new.assay.name = "SCT",      # keep normalized data in SCT assay
  verbose = FALSE
)
message("After SCTransform:")

# saveRDS(Tooth_subset, file = file.path(proc_dir, "Tooth_SCT.rds"))
# Tooth_subset <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Tooth_SCT.rds")

# 3️⃣ Run PCA on SCT assay
Tooth_subset <- RunPCA(
  object  = Tooth_subset,
  assay   = "SCT",
  verbose = FALSE
)
message("After RunPCA:")

# 5️⃣ Run UMAP on selected PCs
Tooth_subset <- RunUMAP(
  object  = Tooth_subset,
  assay   = "SCT",
  dims    = 1:30,
  verbose = FALSE
)

# 6️⃣ Visualize clusters
T_sct <- DimPlot(Tooth_subset, label = TRUE, group.by = "cluster_label", reduction = "umap")

########################################################
# Log Normalisation
######################################################## 
# 2. Log‑normalize the data
Tooth_subset1 <- NormalizeData(
  object = Tooth_subset,
  assay  = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)
message("After NormalizeData:")
print(table(Idents(Tooth_subset1)))       # Show cluster identity counts

# 3. Identify highly variable features
Tooth_subset1 <- FindVariableFeatures(
  object = Tooth_subset1,
  assay  = "RNA",
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)
message("After FindVariableFeatures:")

# 4. Scale data (all genes or HVGs only)
Tooth_subset1 <- ScaleData(
  object = Tooth_subset1,
  assay = "RNA",
  features = VariableFeatures(Tooth_subset1),
  verbose = FALSE
)
message("After ScaleData:")

# 5. Run PCA
Tooth_subset1 <- RunPCA(
  object = Tooth_subset1,
  assay = "RNA",
  features = VariableFeatures(Tooth_subset1),
  verbose = FALSE
)
message("After RunPCA:")

# (From the elbow plot, pick a cutoff; e.g. use dims 1:20 below)

# 7. Run UMAP on selected PCs
Tooth_subset1 <- RunUMAP(
  object = Tooth_subset1,
  assay  = "RNA",
  dims   = 1:30,
  verbose = FALSE
)

DefaultAssay(Tooth_subset1) <- "RNA"       # Set default assay to RNA

T_lognorm <- DimPlot(Tooth_subset1, label = TRUE, group.by = "cluster_label", reduction = "umap")

########################################################
# SubcellType Labelling
########################################################
# 1. Read your annotation table
annotations <- read.csv("/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE184749/Tooth1_Epi_cell_annotations.csv", stringsAsFactors = FALSE)

# 2. Make a named vector: names are the CSV 'cell's, values are the 'new_cluster's
cluster_map <- setNames(annotations$new_cluster, annotations$cell)

# 3. Create (or overwrite) the new_cluster column in meta.data by matching
Tooth_subset1@meta.data$new_cluster <- cluster_map[ Tooth_subset1@meta.data$cell_id ]

# 4️⃣ Find the unmatched cell IDs
unmatched_cells <- setdiff(names(cluster_map), Tooth_subset1@meta.data$cell_id)  # Identify unmatched cell IDs

# 5️⃣ Print or save the list
length(unmatched_cells)        # how many didn’t match
View(unmatched_cells)          # View unmatched cell IDs  

# 📄 Convert to a dataframe
unmatched_df <- data.frame(unmatched_cell_id = unmatched_cells)

# 👀 View in RStudio's data viewer
View(unmatched_df)

# 4. Check it
table( !is.na(Tooth_subset1@meta.data$new_cluster) )  # Count matched cluster labels
head( Tooth_subset1@meta.data[, c("cell_id","new_cluster")] )  # Preview mapping
View(Tooth_subset1@meta.data)                            # View metadata

# grab the two vectors
nc <- Tooth_subset1@meta.data$new_cluster      # new_cluster vector
ct <- Tooth_subset1@meta.data$cluster_label    # original cluster_label vector

# find “empty” entries: NA or empty string
empty_idx <- is.na(nc) | nc == ""

# replace them with the celltype
nc[empty_idx] <- ct[empty_idx]

# put it back
Tooth_subset1@meta.data$new_cluster <- nc

# sanity‐check
table(empty_idx)               # how many got replaced
head( Tooth_subset1@meta.data[, c("cluster_label","new_cluster")] )

View(Tooth_subset1@meta.data)   # Final check of metadata

saveRDS(Tooth_subset, file = file.path(proc_dir, "Tooth_SCT.rds"))       # Save SCTransform processed object
saveRDS(Tooth_subset1, file = file.path(proc_dir, "Tooth_LogNorm.rds"))  # Save LogNormalize processed object

########################################################
# Epithelial Sub-Annotation
########################################################

Tooth_subset <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Tooth_SCT.rds")
# Tooth_subset1 <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Tooth_LogNorm.rds")

T_SCT <- DimPlot(Tooth_subset, label = TRUE, group.by = "cluster_label", reduction = "umap")   # UMAP plot by cluster_label
LEF1 <-FeaturePlot(Tooth_subset, features = c("LEF1"), reduction = "umap", cols = c("lightgrey", "darkred"))  # Feature plot for LEF1 gene

# Subset to only Dental_epithelium cluster
Tooth_epithelial_subset <- subset(
  Tooth_subset,
  subset = cluster_label %in% c("Dental_epithelium")
)

DefaultAssay(Tooth_epithelial_subset)   # Check default assay

# ⃣ Run PCA on SCT assay for epithelial subset
Tooth_epithelial_subset <- RunPCA(
  object  = Tooth_epithelial_subset,
  assay   = "SCT",
  verbose = FALSE
)
message("After RunPCA:")

# Run UMAP on selected PCs for epithelial subset
Tooth_epithelial_subset <- RunUMAP(
  object  = Tooth_epithelial_subset,
  assay   = "SCT",
  dims    = 1:30,
  verbose = FALSE
)

# Visualize epithelial subset UMAP grouped by new_cluster, with nice formatting
DimPlot(
  Tooth_epithelial_subset,
  reduction = "umap",
  group.by = "new_cluster",
  label = TRUE,
  repel = TRUE,                # Prevent label overlap
  label.size = 5
) +
  coord_fixed() +             # Fix 1:1 UMAP axis ratio
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Another UMAP plot grouped by new_cluster
DimPlot(Tooth_epithelial_subset, label = TRUE, group.by = "new_cluster", reduction = "umap")

# Select specific epithelial clusters of interest for subsetting
selected_clusters <- c("EK", "SII", "CL", "OEE", "SIO", "IEE", "SRI", "PA", "SRO", "AM", "SAM")
selected_clusters <- c("EK")  # Override selection to only EK cluster

# Subset to selected epithelial clusters
subset_epithelial <- subset(
  Tooth_epithelial_subset,
  subset = new_cluster %in% selected_clusters
)

# Plot UMAP of subset epithelial clusters with labels and formatting
a <- DimPlot(
  subset_epithelial,
  reduction = "umap",
  group.by = "new_cluster",
  label = TRUE,
  repel = TRUE,                # Prevent label overlap
  label.size = 5
) + 
  coord_fixed() +             # Fix 1:1 UMAP axis ratio
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# Feature plots for marker genes in the EK subset, colored from lightgrey to darkred
EK <-FeaturePlot(subset_epithelial, features = c("SHH", "LEF1", "FGFR2", "CADM1", "PCDH7"), reduction = "umap", cols = c("lightgrey", "darkred"))

# Combine marker feature plots with cluster UMAP plot side by side
EK | a