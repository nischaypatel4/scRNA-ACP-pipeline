#---------------------------------------------
#---------------------------------------------
# Periodontium Dataset QC,Filtering and Integration
# Input : Perio1,2,3,4,5
# Output: Integrated Dataset
#---------------------------------------------
#---------------------------------------------

#---------------------------------------------
# Global seed (Run again, if you repeat the RunUMAP, RunPCA, etc steps)
set.seed(9867) # Set random seed for reproducibility
#---------------------------------------------

#---------------------------------------------
# Load required libraries
#---------------------------------------------
library(Seurat)     # For single-cell analysis
library(dplyr)      # Data manipulation
library(ggplot2)    # Plotting
library(patchwork)  # Combining plots
library(Matrix)     # Sparse matrix operations
library(scales)     # Scale functions for ggplot2
library(DropletUtils) # Droplet-based scRNA-seq utilities
library(harmony)    # Batch correction/integration
library(sctransform) # Normalization and variance stabilization

#---------------------------------------------
# Paths must be set first
#---------------------------------------------
BASE_DIRS <- list( # List of directories for each Perio dataset
  Perio1 = "/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE161266/Perio1",
  Perio2 = "/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE161266/Perio2",
  Perio3 = "/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE161266/Perio3",
  Perio4 = "/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE161266/Perio4",
  Perio5 = "/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE161266/Perio5"
)
proc_dir        <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed"  # Processed data directory
qc_dir          <- "/conglilab/shared/projects/personal_projects/nischay25/project/results/Perio_QC_Filtering_Integration" # QC results directory

dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)  # Create processed dir if not exists
dir.create(qc_dir,   showWarnings = FALSE, recursive = TRUE)  # Create QC results dir if not exists

more_breaks <- function(x) pretty(x, n = 10)  # Helper function for plotting breaks on axes

for (name in names(BASE_DIRS)) {  # Loop through each dataset
  message("\n=== Processing ", name, " ===") # Print current dataset being processed
  
  # 1. Load + create object
  counts <- Read10X(data.dir = BASE_DIRS[[name]])  # Read 10X data
  obj    <- CreateSeuratObject(counts = counts, project = name,
                               min.cells = 3, min.features = 200)  # Create Seurat object with filtering
  
  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-") # Calculate mitochondrial gene percentage
  
  # 2. Violin plots
  v1 <- VlnPlot(obj, "nFeature_RNA", pt.size = .1) +  # Violin plot for number of features
    scale_y_continuous(breaks = more_breaks, labels = comma) +
    ggtitle(paste0(name, ": nFeature_RNA"))
  
  v2 <- VlnPlot(obj, "nCount_RNA",   pt.size = .1) +  # Violin plot for counts
    scale_y_continuous(breaks = more_breaks, labels = comma) +
    ggtitle(paste0(name, ": nCount_RNA"))
  
  v3 <- VlnPlot(obj, "percent.mt",   pt.size = .1) +  # Violin plot for mitochondrial percentage
    scale_y_continuous(breaks = more_breaks, labels = comma) +
    ggtitle(paste0(name, ": percent.mt"))
  
  ggsave(file.path(qc_dir, paste0(name, "_violin.png")), v1 | v2 | v3,
         width = 12, height = 4)  # Save violin plots combined
  
  # 3. Barcode-rank plot
  br    <- barcodeRanks(GetAssayData(obj, slot = "counts"))  # Calculate barcode ranks
  df_br <- data.frame(rank = br$rank, total = br$total)      # Create dataframe for plotting
  p_br  <- ggplot(df_br, aes(rank, total)) +                 # Plot barcode rank curve
    geom_line() +
    geom_hline(yintercept = br@metadata$knee,      linetype = "dashed") +  # Knee threshold
    geom_hline(yintercept = br@metadata$inflection, linetype = "dotted") + # Inflection threshold
    scale_x_log10(labels = comma) + scale_y_log10(labels = comma) +
    ggtitle(paste0(name, ": BarcodeRank vs UMI"))
  
  ggsave(file.path(qc_dir, paste0(name, "_barcode.png")), p_br,
         width = 6, height = 4)  # Save barcode rank plot
  
  # 4. Histograms
  h_u <- ggplot(obj@meta.data, aes(nCount_RNA)) +  # Histogram for UMI counts
    geom_histogram(bins = 50) +
    labs(title = paste0(name, " nUMI")) +
    theme_minimal() +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma)
  
  h_g <- ggplot(obj@meta.data, aes(nFeature_RNA)) +  # Histogram for gene counts
    geom_histogram(bins = 50) +
    labs(title = paste0(name, " nGene")) +
    theme_minimal() +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma)
  
  ggsave(file.path(qc_dir, paste0(name, "_hist.png")), h_u + h_g,
         width = 10, height = 4)  # Save histogram plots side by side
  
  # 5. Scatter plots
  s1 <- FeatureScatter(obj, "nCount_RNA", "percent.mt")  # Scatter plot: UMI counts vs mitochondrial %
  s2 <- FeatureScatter(obj, "nCount_RNA", "nFeature_RNA") # Scatter plot: UMI counts vs gene counts
  
  ggsave(file.path(qc_dir, paste0(name, "_scatter.png")), s1 + s2,
         width = 10, height = 4)  # Save scatter plots side by side
}

# Load each dataset separately (again, maybe for merging)
perio1 <- Read10X(data.dir = "/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE161266/Perio1")
perio1 <- CreateSeuratObject(counts = perio1,   project = "perio1",min.cells = 3,min.features = 200)

perio2 <- Read10X(data.dir = "/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE161266/Perio2")
perio2 <- CreateSeuratObject(counts = perio2,   project = "perio2",min.cells = 3,min.features = 200)

perio3 <- Read10X(data.dir = "/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE161266/Perio3")
perio3 <- CreateSeuratObject(counts = perio3,   project = "perio3",min.cells = 3,min.features = 200)

perio4 <- Read10X(data.dir = "/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE161266/Perio4")
perio4 <- CreateSeuratObject(counts = perio4,   project = "perio4",min.cells = 3,min.features = 200)

perio5 <- Read10X(data.dir = "/conglilab/shared/projects/personal_projects/nischay25/project/data/Ref_Datasets/GSE161266/Perio5")
perio5 <- CreateSeuratObject(counts = perio5,   project = "perio5",min.cells = 3,min.features = 200)

# Merge all datasets into one Seurat object with cell ids preserved
all_perio<-merge(x=perio1, y=c(perio2,perio3,perio4,perio5),add.cell.ids=c("perio1","perio2","perio3","perio4","perio5"))

# store mitochondrial percentage in object meta data
all_perio <- PercentageFeatureSet(all_perio, pattern = "^MT-", col.name = "percent.mt")

# Backup original merged object
all_perio_backup<-all_perio

# Filter cells: keep cells with mitochondrial % < 20 and total counts < 50,000
all_perio<-subset(all_perio, subset = percent.mt < 20 & nCount_RNA < 50000)

# saveRDS(all_perio, file = file.path(proc_dir, "PERIO_Merged_Raw.rds"))  # Optionally save raw merged object

# run sctransform normalization and regression of mitochondrial percentage
set.seed(20191204)  # Set seed for reproducibility
all_perio <- SCTransform(all_perio, vars.to.regress = "percent.mt", verbose = TRUE)

# saveRDS(all_perio,file = "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/PERIO_Merged_SCT.rds") # Optionally save SCT object

# Run PCA on normalized data
all_perio <- RunPCA(all_perio, ndims = 50, verbose = FALSE)

# Run UMAP dimensionality reduction
all_perio <- RunUMAP(all_perio, dims = 1:50, verbose = FALSE)

# Find neighbors for clustering
all_perio <- FindNeighbors(all_perio, dims = 1:50, verbose = FALSE)

# Find clusters
all_perio <- FindClusters(all_perio, verbose = FALSE)

# Plot clusters on UMAP with labels but no legend
DimPlot(all_perio, label = TRUE) + NoLegend()

# Before running Harmony, save pre-Harmony clusters
all_perio$clusters_pre_harmony <- Idents(all_perio)  # Save original clustering

# These lines seem inconsistent but preserve as is
Idents(object = all_perio) <- "clustering_pre_harmony"
Idents(object = all_perio) <- "orig.ident"

# Run Harmony batch correction integration using orig.ident as batch variable
set.seed(20191204)
all_perio<-RunHarmony(all_perio,"orig.ident", theta = 2, plot_convergence = TRUE, nclust = 50, 
                      max.iter.cluster = 20, max.iter.harmony = 5)

# Run UMAP on Harmony embeddings
set.seed(20191204)
all_perio <- RunUMAP(all_perio, dims = 1:50, verbose = FALSE,reduction = "harmony")

# Find neighbors and clusters on Harmony reduction
all_perio <- FindNeighbors(all_perio, dims = 1:50, verbose = FALSE,reduction = "harmony")
all_perio <- FindClusters(all_perio, verbose = FALSE,reduction="harmony")

# Plot Harmony integrated clusters on UMAP
p <- DimPlot(all_perio, label = TRUE,reduction = "umap")
ggsave(file.path(qc_dir, "perio_clustered.png"), p, width = 12, height = 6)

# Feature plots for marker genes to visualize expression patterns
q <- FeaturePlot(
  all_perio,
  features = c("FRZB", "KRT14", "MDK", "FDCSP", "EMCN", "PTPRC"),
  reduction = "umap",
  order = TRUE
)
ggsave(file.path(qc_dir, "perio_markers.png"), q, width = 12, height = 6)

# Custom function to make feature plots with custom colors
make_plot <- function(gene, colors) {
  FeaturePlot(
    object = all_perio,
    features = gene,
    reduction = "umap",
    order = TRUE
  ) +
    scale_color_gradientn(
      colors = colors,
      limits = c(1, 3),
      oob = scales::squish,
      name = NULL
    ) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold.italic", size = 12, hjust = 0.5),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.2, "cm"),
      legend.position = "right",
      plot.margin = margin(2, 2, 2, 2)
    )
}

# Define colors for each gene's feature plot
color_pal <- list(
  FRZB = c("grey90", "pink", "darkred"),
  KRT14 = c("grey90", "gold", "darkorange"),
  MDK = c("grey90", "lightblue", "navy"),
  FDCSP = c("grey90", "wheat", "darkorange"),
  EMCN = c("grey90", "lightgreen", "darkgreen"),
  PTPRC = c("grey90", "plum1", "purple4")
)

# Create individual gene feature plots
p1 <- make_plot("FRZB", color_pal$FRZB)
p2 <- make_plot("KRT14", color_pal$KRT14)
p3 <- make_plot("MDK", color_pal$MDK)
p4 <- make_plot("FDCSP", color_pal$FDCSP)
p5 <- make_plot("EMCN", color_pal$EMCN)
p6 <- make_plot("PTPRC", color_pal$PTPRC)

# Combine all marker plots into 2x3 grid layout
final_plot <- (p1 | p2 ) / (p3 | p4) / (p5 | p6)
ggsave(file.path(qc_dir, "perio_markers.png"), final_plot, width = 12, height = 6)

####PLOT EPITHELIAL MARKER FEATURE PLOTS

###Initial Annotation
# Map cluster IDs to known cell types based on SCT clustering at resolution 0.8
cluster_ids_perio <- as.character(all_perio$SCT_snn_res.0.8)

celltype_annotation_perio <- dplyr::case_when(
  cluster_ids_perio == "0" ~ "Endo1",                # Cluster 0 = Endothelial type 1
  cluster_ids_perio == "9" ~ "Endo2",                # Cluster 9 = Endothelial type 2
  cluster_ids_perio == "12" ~ "SCCs",                # Cluster 12 = SCCs (possibly squamous)
  cluster_ids_perio %in% c("2", "8") ~ "MSCs",       # Clusters 2 & 8 = Mesenchymal stem cells
  cluster_ids_perio %in% c("6", "10") ~ "Fibroblasts", # Clusters 6 & 10 = Fibroblasts
  cluster_ids_perio %in% c("7", "11") ~ "Immune",    # Clusters 7 & 11 = Immune cells
  cluster_ids_perio %in% c("1") ~ "Epithelial1",     # Cluster 1 = Epithelial subtype 1
  cluster_ids_perio %in% c("3", "4", "5") ~ "Epithelial2", # Clusters 3,4,5 = Epithelial subtype 2
  cluster_ids_perio %in% c("13") ~ "Epithelial3",    # Cluster 13 = Epithelial subtype 3
  TRUE ~ "Unknown"                                   # Others = Unknown
)

# Add cell type labels to metadata
all_perio$celltype <- celltype_annotation_perio

# Plot UMAP colored by annotated cell types with labels
a <- DimPlot(all_perio, label = TRUE, group.by= "celltype", reduction = "umap")
ggsave(file.path(qc_dir, "perio_celltype_annotation.png"), a, width = 12, height = 6)

# Save the integrated Seurat object after QC and integration
saveRDS(all_perio, file = file.path(proc_dir, "PERIO_Integrated_1.rds"))
