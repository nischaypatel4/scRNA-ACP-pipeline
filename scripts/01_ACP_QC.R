#---------------------------------------------
#---------------------------------------------
# ACP Single-Cell Dataset QC and Filtering Script
# Inputs  : ACP1, ACP2, ACP3 (10X Genomics data)
# Outputs : QC plots (violin, histograms, scatter, barcode plots) and filtered Seurat objects
#---------------------------------------------
#---------------------------------------------

#---------------------------------------------
# Load required R libraries for analysis
#---------------------------------------------
library(Seurat)        # for single-cell data handling and analysis
library(dplyr)         # for data manipulation
library(ggplot2)       # for visualization
library(patchwork)     # to combine plots
library(Matrix)        # to work with sparse matrices
library(scales)        # for pretty axis labels in plots
library(DropletUtils)  # for barcode rank calculations
library(harmony)       # (future use) for batch correction/integration

#---------------------------------------------
# Define file paths for raw data and outputs
#---------------------------------------------
BASE_DIR_ACP1 <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/raw/acp1"   # ACP1 raw data directory
BASE_DIR_ACP2 <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/raw/acp2"   # ACP2 raw data directory
BASE_DIR_ACP3 <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/raw/acp3"   # ACP3 raw data directory
PROC_DIR      <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed"  # Processed object save location
QC_DIR        <- "/conglilab/shared/projects/personal_projects/nischay25/project/results/ACP_QC_Filtering"  # QC plots output directory

# Create directories if they don’t exist
dir.create(PROC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(QC_DIR,   showWarnings = FALSE, recursive = TRUE)

#---------------------------------------------
# Helper function: create more granular axis breaks for plots
#---------------------------------------------
more_breaks <- function(x) pretty(x, n = 10)

#---------------------------------------------
# 1. Process ACP1 dataset
#---------------------------------------------
message("\n=== Processing acp1 ===")

# Load raw counts from 10X Genomics directory
acp1_data <- Read10X(data.dir = BASE_DIR_ACP1)

# Create Seurat object (min.cells=3, min.features=200)
acp1_obj  <- CreateSeuratObject(counts = acp1_data, project = "acp1", min.cells = 3, min.features = 200)

# Calculate percentage of mitochondrial genes for each cell
acp1_obj$percent.mt <- PercentageFeatureSet(acp1_obj, pattern = "^MT-")

#---------------------------------------------
# Generate QC plots for ACP1
#---------------------------------------------

# Violin plots for detected genes, UMI counts, and mitochondrial percentage
v1_acp1 <- VlnPlot(acp1_obj, features = "nFeature_RNA", pt.size = 0.1) +
  scale_y_continuous(breaks = more_breaks, labels = comma) + labs(title = "acp1: nFeature_RNA")
v2_acp1 <- VlnPlot(acp1_obj, features = "nCount_RNA", pt.size = 0.1) +
  scale_y_continuous(breaks = more_breaks, labels = comma) + labs(title = "acp1: nCount_RNA")
v3_acp1 <- VlnPlot(acp1_obj, features = "percent.mt", pt.size = 0.1) +
  scale_y_continuous(breaks = more_breaks, labels = comma) + labs(title = "acp1: percent.mt")
gi_acp1 <- (v1_acp1 | v2_acp1 | v3_acp1)
ggsave(file.path(QC_DIR, "acp1_violin.png"), gi_acp1, width = 12, height = 4)

# Barcode rank plot to identify knee/inflection point
br_acp1   <- barcodeRanks(GetAssayData(acp1_obj, slot = "counts"))
df_br1    <- data.frame(rank = br_acp1$rank, total = br_acp1$total)
p_acp1    <- ggplot(df_br1, aes(rank, total)) +
  geom_line() + geom_hline(yintercept = br_acp1@metadata$knee, linetype = "dashed") +
  geom_hline(yintercept = br_acp1@metadata$inflection, linetype = "dotted") +
  scale_x_log10(labels = comma) + scale_y_log10(labels = comma) +
  labs(title = "acp1: BarcodeRank vs UMI")
ggsave(file.path(QC_DIR, "acp1_barcode.png"), p_acp1, width = 6, height = 4)

# Histogram for UMI counts per cell
hist1_u <- ggplot(acp1_obj@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "black") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(title = "acp1 nUMI", x = "nUMI", y = "Cell count") +
  theme_minimal()

# Histogram for number of genes per cell
hist1_g <- ggplot(acp1_obj@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 50, fill = "darkorange", color = "black") +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(title = "acp1 nGene", x = "nGene", y = "Cell count") +
  theme_minimal()
hist1   <- hist1_u + hist1_g
ggsave(file.path(QC_DIR, "acp1_hist.png"), hist1, width = 10, height = 4)

# Scatter plots: nUMI vs percent.mt and nUMI vs nGene
sc1_1 <- FeatureScatter(acp1_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
sc1_2 <- FeatureScatter(acp1_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
scat1 <- sc1_1 + sc1_2
ggsave(file.path(QC_DIR, "acp1_scatter.png"), scat1, width = 10, height = 4)

#---------------------------------------------
# Filter ACP1 dataset based on QC thresholds
#---------------------------------------------

# Define cell-level QC cutoffs
md1       <- acp1_obj@meta.data
keep1     <- md1$nFeature_RNA > 1200 & md1$nCount_RNA < 45000 & md1$percent.mt < 25

# Subset cells passing QC
acp1_filt <- subset(acp1_obj, cells = rownames(md1)[keep1])

# Filter genes: retain only genes expressed in ≥20 cells
g1_counts <- GetAssayData(acp1_filt, slot = "counts")
g1_frac   <- Matrix::rowSums(g1_counts > 0)
keep1_g   <- names(g1_frac[g1_frac >= 20])

# Create final filtered object
acp1_final<- subset(acp1_filt, features = keep1_g)

#---------------------------------------------
# Summarize ACP1 filtering results
#---------------------------------------------

# Cell and gene counts (raw, filtered, final)
cat(sprintf("acp1_raw: %d cells × %d genes\n", ncol(acp1_obj), nrow(acp1_obj)))
cat(sprintf("acp1_filtered: %d cells × %d genes\n", ncol(acp1_filt), nrow(acp1_filt)))
cat(sprintf("acp1_final: %d cells × %d genes\n", ncol(acp1_final), nrow(acp1_final)))

# QC metric ranges before filtering
umi_raw   <- acp1_obj@meta.data$nCount_RNA
gene_raw  <- acp1_obj@meta.data$nFeature_RNA
mito_raw  <- acp1_obj@meta.data$percent.mt
cat(sprintf("acp1_raw nUMI range: %.0f - %.0f\n", min(umi_raw), max(umi_raw)))
cat(sprintf("acp1_raw nGene range: %.0f - %.0f\n", min(gene_raw), max(gene_raw)))
cat(sprintf("acp1_raw percent.mt range: %.2f%% - %.2f%%\n", min(mito_raw), max(mito_raw)))

# QC metric ranges after filtering
umi_final  <- acp1_final@meta.data$nCount_RNA
gene_final <- acp1_final@meta.data$nFeature_RNA
mito_final <- acp1_final@meta.data$percent.mt
cat(sprintf("acp1_final nUMI range: %.0f - %.0f\n", min(umi_final), max(umi_final)))
cat(sprintf("acp1_final nGene range: %.0f - %.0f\n", min(gene_final), max(gene_final)))
cat(sprintf("acp1_final percent.mt range: %.2f%% - %.2f%%\n", min(mito_final), max(mito_final)))

# Save ACP1 filtered objects
saveRDS(acp1_filt,  file.path(PROC_DIR, "acp1_filtered.rds"))
saveRDS(acp1_final,file.path(PROC_DIR, "acp1_final.rds"))

#---------------------------------------------
# 2. Process ACP2 dataset
#---------------------------------------------
message("\n=== Processing acp2 ===")

# Load ACP2 counts and create Seurat object
acp2_data <- Read10X(data.dir = BASE_DIR_ACP2)
acp2_obj  <- CreateSeuratObject(counts = acp2_data, project = "acp2", min.cells = 3, min.features = 200)
acp2_obj$percent.mt <- PercentageFeatureSet(acp2_obj, pattern = "^MT-")

# Violin plots for ACP2
v1_acp2 <- VlnPlot(acp2_obj, features = "nFeature_RNA", pt.size = 0.1) + scale_y_continuous(breaks = more_breaks, labels = comma) + labs(title = "acp2: nFeature_RNA")
v2_acp2 <- VlnPlot(acp2_obj, features = "nCount_RNA", pt.size = 0.1) + scale_y_continuous(breaks = more_breaks, labels = comma) + labs(title = "acp2: nCount_RNA")
v3_acp2 <- VlnPlot(acp2_obj, features = "percent.mt", pt.size = 0.1) + scale_y_continuous(breaks = more_breaks, labels = comma) + labs(title = "acp2: percent.mt")
gi_acp2 <- (v1_acp2 | v2_acp2 | v3_acp2)
ggsave(file.path(QC_DIR, "acp2_violin.png"), gi_acp2, width = 12, height = 4)

# Barcode rank plot ACP2
br_acp2<- barcodeRanks(GetAssayData(acp2_obj, slot = "counts"))
df_br2 <- data.frame(rank = br_acp2$rank, total = br_acp2$total)
p_acp2 <- ggplot(df_br2, aes(rank, total)) + geom_line() + geom_hline(yintercept = br_acp2@metadata$knee, linetype = "dashed") + geom_hline(yintercept = br_acp2@metadata$inflection, linetype = "dotted") + scale_x_log10(labels = comma) + scale_y_log10(labels = comma) + labs(title = "acp2: BarcodeRank vs UMI")
ggsave(file.path(QC_DIR, "acp2_barcode.png"), p_acp2, width = 6, height = 4)

# Histograms for ACP2
hist2_u <- ggplot(acp2_obj@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 500, fill = "steelblue", color = "black") +
  scale_x_continuous(breaks = more_breaks) +
  scale_y_continuous(breaks = more_breaks) +
  labs(title = "acp2: Distribution of nUMI", x = "nUMI", y = "Cell Count") +
  theme_minimal()

hist2_g <- ggplot(acp2_obj@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill = "darkorange", color = "black") +
  scale_x_continuous(breaks = more_breaks) +
  scale_y_continuous(breaks = more_breaks) +
  labs(title = "acp2: Distribution of nGene", x = "nGene", y = "Cell Count") +
  theme_minimal()

hist2 <- hist2_u + hist2_g

# Scatter plots for ACP2
sc2_1 <- FeatureScatter(acp2_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
sc2_2 <- FeatureScatter(acp2_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
scat2  <- sc2_1 + sc2_2
ggsave(file.path(QC_DIR, "acp2_scatter.png"), scat2, width = 10, height = 4)

# Filter ACP2 cells based on QC thresholds
md2        <- acp2_obj@meta.data
keep2      <- md2$nFeature_RNA > 1700 & md2$nCount_RNA < 45000 & md2$percent.mt < 25
acp2_filt  <- subset(acp2_obj, cells = rownames(md2)[keep2])

# Filter ACP2 genes (≥20 cells)
g2_counts  <- GetAssayData(acp2_filt, slot = "counts")
g2_frac    <- Matrix::rowSums(g2_counts > 0)
keep2_g    <- names(g2_frac[g2_frac >= 20])
acp2_final <- subset(acp2_filt, features = keep2_g)

# Summary ACP2
cat(sprintf("acp2_raw: %d cells × %d genes\n", ncol(acp2_obj), nrow(acp2_obj)))
cat(sprintf("acp2_filtered: %d cells × %d genes\n", ncol(acp2_filt), nrow(acp2_filt)))
cat(sprintf("acp2_final: %d cells × %d genes\n", ncol(acp2_final), nrow(acp2_final)))

# QC metric ranges (ACP2 raw)
umi_raw   <- acp2_obj@meta.data$nCount_RNA
gene_raw  <- acp2_obj@meta.data$nFeature_RNA
mito_raw  <- acp2_obj@meta.data$percent.mt
cat(sprintf("acp2_raw nUMI range: %.0f - %.0f\n", min(umi_raw), max(umi_raw)))
cat(sprintf("acp2_raw nGene range: %.0f - %.0f\n", min(gene_raw), max(gene_raw)))
cat(sprintf("acp2_raw percent.mt range: %.2f%% - %.2f%%\n", min(mito_raw), max(mito_raw)))

# QC metric ranges (ACP2 filtered)
umi_final  <- acp2_final@meta.data$nCount_RNA
gene_final <- acp2_final@meta.data$nFeature_RNA
mito_final <- acp2_final@meta.data$percent.mt
cat(sprintf("acp2_final nUMI range: %.0f - %.0f\n", min(umi_final), max(umi_final)))
cat(sprintf("acp2_final nGene range: %.0f - %.0f\n", min(gene_final), max(gene_final)))
cat(sprintf("acp2_final percent.mt range: %.2f%% - %.2f%%\n", min(mito_final), max(mito_final)))

# Save ACP2 filtered objects
saveRDS(acp2_filt,  file.path(PROC_DIR, "acp2_filtered.rds"))
saveRDS(acp2_final,file.path(PROC_DIR, "acp2_final.rds"))

#---------------------------------------------
# 3. Process ACP3 dataset
#---------------------------------------------
message("\n=== Processing acp3 ===")

# Load ACP3 counts and create Seurat object
acp3_data <- Read10X(data.dir = BASE_DIR_ACP3)
acp3_obj  <- CreateSeuratObject(counts = acp3_data, project = "acp3", min.cells = 3, min.features = 200)
acp3_obj$percent.mt <- PercentageFeatureSet(acp3_obj, pattern = "^MT-")

# Violin plots ACP3
v1_acp3 <- VlnPlot(acp3_obj, features = "nFeature_RNA", pt.size = 0.1) + scale_y_continuous(breaks = more_breaks, labels = comma) + labs(title = "acp3: nFeature_RNA")
v2_acp3 <- VlnPlot(acp3_obj, features = "nCount_RNA", pt.size = 0.1) + scale_y_continuous(breaks = more_breaks, labels = comma) + labs(title = "acp3: nCount_RNA")
v3_acp3 <- VlnPlot(acp3_obj, features = "percent.mt", pt.size = 0.1) + scale_y_continuous(breaks = more_breaks, labels = comma) + labs(title = "acp3: percent.mt")
gi_acp3 <- (v1_acp3 | v2_acp3 | v3_acp3)
ggsave(file.path(QC_DIR, "acp3_violin.png"), gi_acp3, width = 12, height = 4)

# Barcode rank plot ACP3
br_acp3<- barcodeRanks(GetAssayData(acp3_obj, slot = "counts"))
df_br3 <- data.frame(rank = br_acp3$rank, total = br_acp3$total)
p_acp3 <- ggplot(df_br3, aes(rank, total)) +
  geom_line() +
  geom_hline(yintercept = br_acp3@metadata$knee, linetype = "dashed") +
  geom_hline(yintercept = br_acp3@metadata$inflection, linetype = "dotted") +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  labs(title = "acp3: BarcodeRank vs UMI")
ggsave(file.path(QC_DIR, "acp3_barcode.png"), p_acp3, width = 6, height = 4)

# Histograms for ACP3: UMI counts and number of genes
hist3_u <- ggplot(acp3_obj@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(binwidth = 500, fill = "steelblue", color = "black") +
  scale_x_continuous(breaks = more_breaks) +
  scale_y_continuous(breaks = more_breaks) +
  labs(title = "acp3: Distribution of nUMI", x = "nUMI", y = "Cell Count") +
  theme_minimal()

hist3_g <- ggplot(acp3_obj@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(binwidth = 100, fill = "darkorange", color = "black") +
  scale_x_continuous(breaks = more_breaks) +
  scale_y_continuous(breaks = more_breaks) +
  labs(title = "acp3: Distribution of nGene", x = "nGene", y = "Cell Count") +
  theme_minimal()

hist3 <- hist3_u + hist3_g
ggsave(file.path(QC_DIR, "acp3_hist.png"), hist3, width = 10, height = 4)

# Scatter plots for ACP3 QC metrics
sc3_1 <- FeatureScatter(acp3_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
sc3_2 <- FeatureScatter(acp3_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
scat3  <- sc3_1 + sc3_2
ggsave(file.path(QC_DIR, "acp3_scatter.png"), scat3, width = 10, height = 4)

#---------------------------------------------
# Filter ACP3 dataset based on QC thresholds
#---------------------------------------------

# Cell-level QC filter for ACP3
md3        <- acp3_obj@meta.data
keep3      <- md3$nFeature_RNA > 1100 & md3$nCount_RNA < 50000 & md3$percent.mt < 25
acp3_filt  <- subset(acp3_obj, cells = rownames(md3)[keep3])

# Gene filter for ACP3: keep genes expressed in ≥20 cells
g3_counts  <- GetAssayData(acp3_filt, slot = "counts")
g3_frac    <- Matrix::rowSums(g3_counts > 0)
keep3_g    <- names(g3_frac[g3_frac >= 20])
acp3_final <- subset(acp3_filt, features = keep3_g)

#---------------------------------------------
# Summarize ACP3 filtering results
#---------------------------------------------

# Cell and gene counts for ACP3
cat(sprintf("acp3_raw: %d cells × %d genes\n", ncol(acp3_obj), nrow(acp3_obj)))
cat(sprintf("acp3_filtered: %d cells × %d genes\n", ncol(acp3_filt), nrow(acp2_filt)))   # (NOTE: probably a typo using acp2_filt here)
cat(sprintf("acp3_final: %d cells × %d genes\n", ncol(acp3_final), nrow(acp3_final)))

# QC metric ranges before filtering (raw ACP3)
umi_raw   <- acp3_obj@meta.data$nCount_RNA
gene_raw  <- acp3_obj@meta.data$nFeature_RNA
mito_raw  <- acp3_obj@meta.data$percent.mt
cat(sprintf("acp3_raw nUMI range: %.0f - %.0f\n", min(umi_raw), max(umi_raw)))
cat(sprintf("acp3_raw nGene range: %.0f - %.0f\n", min(gene_raw), max(gene_raw)))
cat(sprintf("acp3_raw percent.mt range: %.2f%% - %.2f%%\n", min(mito_raw), max(mito_raw)))

# QC metric ranges after filtering (ACP3)
umi_final  <- acp3_final@meta.data$nCount_RNA
gene_final <- acp3_final@meta.data$nFeature_RNA
mito_final <- acp3_final@meta.data$percent.mt
cat(sprintf("acp3_final nUMI range: %.0f - %.0f\n", min(umi_final), max(umi_final)))
cat(sprintf("acp3_final nGene range: %.0f - %.0f\n", min(gene_final), max(gene_final)))
cat(sprintf("acp3_final percent.mt range: %.2f%% - %.2f%%\n", min(mito_final), max(mito_final)))

# Save ACP3 filtered Seurat objects
saveRDS(acp3_filt,  file.path(PROC_DIR, "acp3_filtered.rds"))
saveRDS(acp3_final,file.path(PROC_DIR, "acp3_final.rds"))

#---------------------------------------------
# Final message after processing all samples
#---------------------------------------------
message("\nAll samples processed successfully.")

#---------------------------------------------
# Merge ACP1, ACP2, ACP3 after QC and filtering
# and generate a final violin plot
#---------------------------------------------

combined <- merge(
  x = acp1_final,
  y = list(acp2_final, acp3_final),
  add.cell.ids = c("acp1", "acp2", "acp3"),
  project = "acp_combined"
)

# Violin plot of key QC metrics across all filtered datasets
ACP_violin <- VlnPlot(combined, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"))
ggsave(file.path(QC_DIR, "acp_violin_afterQC.png"), ACP_violin , width = 10, height = 4)