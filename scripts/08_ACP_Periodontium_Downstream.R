#---------------------------------------------
#---------------------------------------------
# ACP_Periondium Downstream Analysis
# Input : Integrated Objects
# Output: Gene Ontology
#---------------------------------------------
#---------------------------------------------

#-------------------------------------------------------------------------------
# Libraries and Loading
#-------------------------------------------------------------------------------

library(Seurat)         # For single-cell data processing
library(dplyr)          # Data manipulation
library(ggplot2)        # Plotting
library(patchwork)      # Combining plots
library(Matrix)         # Matrix operations
library(scales)         # Scales for plotting
library(DropletUtils)   # Droplet-based scRNA-seq tools
library(harmony)        # Batch correction and integration
library(RColorBrewer)   # Color palettes
library(tibble)         # Enhanced data frames
library(clusterProfiler) # Functional enrichment analysis
library(ReactomePA)     # Reactome pathway analysis
library(org.Hs.eg.db)   # Human gene annotation database
library(enrichplot)     # Enrichment result visualization

proc_dir <- "/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Downstream"  # Output directory

#-------------------------------------------------------------------------------
# Data Loading
#-------------------------------------------------------------------------------
# Load Seurat objects containing raw expression data
acp_Epi_KE_raw <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Downstream/acp_Epi_KE_raw.rds")
acp_Epi_PE_raw <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Downstream/acp_Epi_PE_raw.rds")

perio_Epi1_raw <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Downstream/perio_Epi1_raw.rds")
perio_Epi2_raw <- readRDS("/conglilab/shared/projects/personal_projects/nischay25/project/data/processed/Downstream/perio_Epi2_raw.rds")

#-------------------------------------------------------------------------------
# Normalisation
#-------------------------------------------------------------------------------
# 1. Merge ACP‐raw objects for KE vs PE ----
acl_KE_PE <- merge(
  x = acp_Epi_KE_raw,
  y = acp_Epi_PE_raw,
  add.cell.ids = c("KE","PE"),       # Add cell ID prefixes for source
  project = "ACP_KE_PE"              # Project name
)

# 2. Normalize the merged ACP object ----
DefaultAssay(acl_KE_PE) <- "RNA"    # Set assay to RNA for normalization
acl_KE_PE <- NormalizeData(
  object = acl_KE_PE,
  normalization.method = "LogNormalize",  # Log-normalize counts
  scale.factor = 1e4,
  verbose = FALSE
)

# 3. Run DE for KE vs PE ----
# Make sure your metadata has the correct cell‐type labels
# (e.g. maybe stored in acl_KE_PE$celltype)
Idents(acl_KE_PE) <- acl_KE_PE$celltype   # Set cluster identities to celltype column

# 1. Join layers if needed (after merge)
acl_KE_PE <- JoinLayers(acl_KE_PE)         # Combine layers if present (e.g., raw + normalized)

# 2. Normalize (if not already done)
DefaultAssay(acl_KE_PE) <- "RNA"
acl_KE_PE <- NormalizeData(acl_KE_PE, normalization.method = "LogNormalize")

# 3. Set identity and run DE
Idents(acl_KE_PE) <- acl_KE_PE$celltype    # Confirm identity class

deg_KE_vs_PE <- FindMarkers(
  object    = acl_KE_PE,
  ident.1   = "acp_Epi_KE",                 # Group 1: KE cells
  ident.2   = "acp_Epi_PE",                 # Group 2: PE cells
  assay     = "RNA",
  slot      = "data",                       # Use normalized data
  test.use  = "wilcox",                     # Wilcoxon test for DE
  min.pct   = 0.1,                         # Minimum fraction of cells expressing gene
  logfc.threshold = 0.25                    # Minimum log fold change threshold
)

# 4. Add gene column and filter
deg_KE_vs_PE$gene <- rownames(deg_KE_vs_PE)   # Add gene names as a column
deg_KE_vs_PE <- deg_KE_vs_PE[deg_KE_vs_PE$p_val_adj < 0.05, ]  # Filter significant DEGs

# 4. Merge Perio‐raw objects for Epi1 vs Epi2 ----
perio_1_2 <- merge(
  x = perio_Epi1_raw,
  y = perio_Epi2_raw,
  add.cell.ids = c("Epi1","Epi2"),          # Add cell ID prefixes for source
  project = "PERIO_Epi1_Epi2"
)

# Join layered RNA data before normalization
perio_1_2 <- JoinLayers(perio_1_2)           # Join layers (raw, normalized)

# Set assay and normalize
DefaultAssay(perio_1_2) <- "RNA"
perio_1_2 <- NormalizeData(
  object = perio_1_2,
  normalization.method = "LogNormalize",
  scale.factor = 1e4,
  verbose = FALSE
)

# 6. Run DE for Perio Epi1 vs Epi2 ----
Idents(perio_1_2) <- perio_1_2$celltype       # Set identities

# Run DE
deg_Epi1_vs_Epi2 <- FindMarkers(
  object    = perio_1_2,
  ident.1   = "perio_Epi1",                    # Group 1: Epi1 cells
  ident.2   = "perio_Epi2",                    # Group 2: Epi2 cells
  assay     = "RNA",
  slot      = "data",
  test.use  = "wilcox",
  min.pct   = 0.1,
  logfc.threshold = 0.25
)

# Add gene column and filter for significance
deg_Epi1_vs_Epi2$gene <- rownames(deg_Epi1_vs_Epi2)   # Add gene names as column
deg_Epi1_vs_Epi2 <- deg_Epi1_vs_Epi2[deg_Epi1_vs_Epi2$p_val_adj < 0.05, ]   # Filter significant

# 7. Intersect the two DEG lists ----
common_genes <- intersect(deg_KE_vs_PE$gene, deg_Epi1_vs_Epi2$gene)   # Find overlapping DE genes
length(common_genes)  # how many overlap?
common_genes         # Print overlapping genes

#-------------------------------------------------------------------------------
# GeneEnrichment
#-------------------------------------------------------------------------------
# 8. GO enrichment on the common genes ----
ego_common <- enrichGO(
  gene          = common_genes,             # Genes for enrichment
  OrgDb         = org.Hs.eg.db,             # Human annotation DB
  keyType       = "SYMBOL",                  # Input gene IDs are SYMBOLs
  ont           = "BP",                      # Ontology: Biological Process
  pAdjustMethod = "BH",                      # Multiple testing correction method
  pvalueCutoff  = 0.05,                      # P-value cutoff
  qvalueCutoff  = 0.1                        # Q-value cutoff
)

# View and plot
head(as.data.frame(ego_common))              # Show top enriched GO terms
dotplot(ego_common, showCategory = 20) +    # Dotplot of top 20 GO terms
  ggtitle("GO BP Enrichment of Common DEGs")
# library(msigdbr)                          # (Commented out, optional)

View(as.data.frame(ego_common))              # View enrichment results in spreadsheet viewer

#-------------------------------------------------------------------------------
# KEGG Analysis
#-------------------------------------------------------------------------------
# Convert SYMBOLs to ENTREZ IDs
gene_df <- bitr(
  common_genes,
  fromType = "SYMBOL",                      # From gene symbol
  toType   = "ENTREZID",                    # To ENTREZ ID
  OrgDb    = org.Hs.eg.db
)
entrez_ids <- gene_df$ENTREZID              # Extract ENTREZ IDs

# 2. Run KEGG enrichment
ekegg <- enrichKEGG(
  gene         = entrez_ids,                 # Gene list
  organism     = 'hsa',                       # Human KEGG organism code
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)

# 3. Inspect results
head(as.data.frame(ekegg))                   # Show top KEGG pathways

# 4. Visualize
dotplot(ekegg, showCategory = 15) +          # Dotplot top 15 KEGG pathways
  ggtitle("KEGG Pathway Enrichment of Common DEGs")


#-------------------------------------------------------------------------------
# ReactomePA
#-------------------------------------------------------------------------------
# Use the same ENTREZ ID list:
rpa <- enrichPathway(
  gene         = entrez_ids,                  # Gene list
  organism     = "human",                      # Organism name
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  readable     = TRUE           # converts back to gene SYMBOLs in the output
)

# View top Reactome pathways
head(as.data.frame(rpa))

# Plot
dotplot(rpa, showCategory = 15) +             # Dotplot top 15 Reactome pathways
  ggtitle("Reactome Pathway Enrichment of Common DEGs")

#-------------------------------------------------------------------------------
# Gene-Concept Network Plot
#-------------------------------------------------------------------------------
# 8.1.1 Gene‑Concept Network (cnetplot)
cnetplot(
  ego_common,
  showCategory = 10,       # Show top 10 GO terms
  foldChange   = NULL,     # Optionally provide avg log₂FC named vector here
  circular     = FALSE,
  colorEdge    = TRUE
) + ggtitle("Gene–GO Term Network")

# 8.1.2 Enrichment Map (emapplot)
emapplot(
  pairwise_termsim(ego_common),     # Calculate similarity between GO terms
  showCategory = 20
) + ggtitle("GO Term Similarity Network")

#-------------------------------------------------------------------------------
# geneCluster GO comparision
#-------------------------------------------------------------------------------
# prepare a named list of gene vectors
gene_list <- list(
  KE_vs_PE    = deg_KE_vs_PE$gene,          # DEG genes from KE vs PE comparison
  Epi1_vs_Epi2 = deg_Epi1_vs_Epi2$gene     # DEG genes from Epi1 vs Epi2 comparison
)

# run compareCluster
cc <- compareCluster(
  geneCluster   = gene_list,
  fun           = "enrichGO",                # Use GO enrichment function
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",                       # Biological process ontology
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# visualize
dotplot(cc, showCategory = 15) +               # Dotplot of GO terms for comparison
  ggtitle("GO BP Enrichment: KE vs PE vs Epi1 vs Epi2")