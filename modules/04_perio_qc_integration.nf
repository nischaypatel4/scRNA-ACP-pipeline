/*
================================================================================
  Module: PERIO_QC
  Purpose : QC, filtering, and Harmony integration of periodontium reference
            dataset (GSE161266: 5 samples)
  Input   : Periodontium raw data directory (GSE161266)
  Output  : Merged raw and Harmony-integrated Seurat objects
================================================================================
*/

process PERIO_QC {
    publishDir "${params.outdir}/Perio_QC_Filtering_Integration", mode: 'copy'

    input:
    path perio_dir

    output:
    path "PERIO_Merged_Raw.rds",  emit: merged_raw
    path "PERIO_Integrated.rds",  emit: integrated
    path "plots/*.png"

    script:
    """
    mkdir -p plots

    Rscript - <<'REOF'
    library(Seurat)
    library(harmony)
    library(ggplot2)
    library(patchwork)
    library(Matrix)
    library(DropletUtils)
    library(dplyr)

    set.seed(${params.seed})

    perio_dir <- "${perio_dir}"
    samples   <- c("Perio1", "Perio2", "Perio3", "Perio4", "Perio5")

    # ── Load and QC each Periodontium sample ──────────────────────────────────
    seurat_list <- lapply(samples, function(s) {
        path <- file.path(perio_dir, s)
        counts <- Read10X(data.dir = path)
        obj <- CreateSeuratObject(counts, project = s,
                                  min.cells = 3, min.features = 200)
        obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
        obj\$sample <- s

        # Filter
        obj <- subset(obj,
                      subset = nFeature_RNA > 200 &
                               nFeature_RNA < 6000 &
                               percent.mt < 20)
        obj
    })
    names(seurat_list) <- samples

    # ── Merge ─────────────────────────────────────────────────────────────────
    merged_raw <- merge(seurat_list[[1]],
                        y = seurat_list[2:5],
                        add.cell.ids = samples)
    saveRDS(merged_raw, "PERIO_Merged_Raw.rds")

    # ── SCTransform + PCA + Harmony ───────────────────────────────────────────
    merged_sct <- SCTransform(merged_raw, vars.to.regress = "percent.mt",
                              verbose = FALSE)
    merged_sct <- RunPCA(merged_sct, npcs = 30, verbose = FALSE)
    integrated <- RunHarmony(merged_sct, group.by.vars = "sample",
                             reduction = "pca", reduction.save = "harmony",
                             verbose = FALSE)
    integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:20,
                          reduction.name = "umap.harmony")
    integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:20)
    integrated <- FindClusters(integrated, resolution = 0.4, algorithm = 4)

    p_umap <- DimPlot(integrated, reduction = "umap.harmony",
                      group.by = c("sample", "seurat_clusters"), ncol = 2)
    ggsave("plots/PERIO_integrated_UMAP.png", p_umap, width = 14, height = 6)

    saveRDS(integrated, "PERIO_Integrated.rds")
    message("Periodontium integration complete: ", ncol(integrated), " cells")
    REOF
    """
}
