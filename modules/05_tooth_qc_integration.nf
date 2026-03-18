/*
================================================================================
  Module: TOOTH_QC
  Purpose : Load, QC, and process the tooth developmental reference dataset
            (GSE184749) using LogNormalize and SCTransform
  Input   : Tooth dataset directory (GSE184749)
  Output  : Processed Seurat objects (LogNorm and SCT)
================================================================================
*/

process TOOTH_QC {
    publishDir "${params.outdir}/Tooth_QC", mode: 'copy'

    input:
    path tooth_dir

    output:
    path "Tooth.rds",         emit: processed
    path "Tooth_LogNorm.rds", emit: lognorm
    path "Tooth_SCT.rds",     emit: sct
    path "plots/*.png"

    script:
    """
    mkdir -p plots

    Rscript - <<'REOF'
    library(Seurat)
    library(Matrix)
    library(ggplot2)
    library(patchwork)
    library(harmony)

    set.seed(${params.seed})

    tooth_dir <- "${tooth_dir}"

    # ── Load count matrix + annotations ──────────────────────────────────────
    count_matrix <- Matrix::readMM(
        file.path(tooth_dir, "GSE184749_all_data_counts.mtx.gz"))
    gene_anno <- read.csv(
        file.path(tooth_dir, "GSE184749_all_data_gene_annotations.csv.gz"))
    cell_anno <- read.csv(
        file.path(tooth_dir, "GSE184749_all_data_cell_annotations.csv.gz"))

    rownames(gene_anno) <- gene_anno[, 1]
    rownames(cell_anno) <- cell_anno[, 1]
    rownames(count_matrix) <- rownames(gene_anno)
    colnames(count_matrix) <- rownames(cell_anno)

    # ── Create Seurat object ──────────────────────────────────────────────────
    tooth <- CreateSeuratObject(counts = count_matrix,
                                project = "Tooth",
                                meta.data = cell_anno)
    tooth[["percent.mt"]] <- PercentageFeatureSet(tooth, pattern = "^MT-")

    # Filter
    tooth <- subset(tooth,
                    subset = nFeature_RNA > 200 &
                             nFeature_RNA < 8000 &
                             percent.mt < 25)

    saveRDS(tooth, "Tooth.rds")

    # ── LogNormalize pipeline ─────────────────────────────────────────────────
    tooth_lognorm <- NormalizeData(tooth)
    tooth_lognorm <- FindVariableFeatures(tooth_lognorm)
    tooth_lognorm <- ScaleData(tooth_lognorm)
    tooth_lognorm <- RunPCA(tooth_lognorm, npcs = 30, verbose = FALSE)
    tooth_lognorm <- RunUMAP(tooth_lognorm, dims = 1:20)
    tooth_lognorm <- FindNeighbors(tooth_lognorm, dims = 1:20)
    tooth_lognorm <- FindClusters(tooth_lognorm, resolution = 0.4)
    saveRDS(tooth_lognorm, "Tooth_LogNorm.rds")

    # ── SCTransform pipeline ──────────────────────────────────────────────────
    tooth_sct <- SCTransform(tooth, vars.to.regress = "percent.mt",
                             verbose = FALSE)
    tooth_sct <- RunPCA(tooth_sct, npcs = 30, verbose = FALSE)
    tooth_sct <- RunUMAP(tooth_sct, dims = 1:20)
    tooth_sct <- FindNeighbors(tooth_sct, dims = 1:20)
    tooth_sct <- FindClusters(tooth_sct, resolution = 0.4)

    p_umap <- DimPlot(tooth_sct, group.by = "seurat_clusters")
    ggsave("plots/Tooth_SCT_UMAP.png", p_umap, width = 8, height = 6)

    saveRDS(tooth_sct, "Tooth_SCT.rds")
    message("Tooth QC complete: ", ncol(tooth_sct), " cells")
    REOF
    """
}
