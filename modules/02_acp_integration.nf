/*
================================================================================
  Module: ACP_INTEGRATION
  Purpose : Merge and integrate ACP1, ACP2, ACP3 using Harmony and SCTransform
  Input   : Three filtered ACP Seurat objects
  Output  : Merged raw, merged SCT, and Harmony-integrated Seurat objects
================================================================================
*/

process ACP_INTEGRATION {
    publishDir "${params.outdir}/ACP_Integration", mode: 'copy'

    input:
    path acp1_rds
    path acp2_rds
    path acp3_rds

    output:
    path "ACP_Merged_Raw.rds",  emit: merged_raw
    path "ACP_Merged_SCT.rds",  emit: merged_sct
    path "ACP_Integrated.rds",  emit: integrated
    path "acp1_final.rds",      emit: acp_final_rds
    path "plots/*.png"

    script:
    """
    mkdir -p plots

    Rscript - <<'REOF'
    library(Seurat)
    library(harmony)
    library(ggplot2)
    library(patchwork)
    library(sctransform)

    set.seed(${params.seed})

    # ── Load filtered objects ─────────────────────────────────────────────────
    acp1 <- readRDS("${acp1_rds}")
    acp2 <- readRDS("${acp2_rds}")
    acp3 <- readRDS("${acp3_rds}")

    # Add sample identifiers
    acp1\$sample <- "ACP1"
    acp2\$sample <- "ACP2"
    acp3\$sample <- "ACP3"

    # ── Merge raw objects ─────────────────────────────────────────────────────
    merged_raw <- merge(acp1, y = list(acp2, acp3),
                        add.cell.ids = c("ACP1", "ACP2", "ACP3"))
    saveRDS(merged_raw, "ACP_Merged_Raw.rds")

    # Copy acp1_final for Scimilarity module
    file.copy("${acp1_rds}", "acp1_final.rds")

    # ── SCTransform normalisation ─────────────────────────────────────────────
    merged_sct <- SCTransform(merged_raw, vars.to.regress = "percent.mt",
                              verbose = FALSE)
    saveRDS(merged_sct, "ACP_Merged_SCT.rds")

    # ── PCA ───────────────────────────────────────────────────────────────────
    merged_sct <- RunPCA(merged_sct, npcs = 30, verbose = FALSE)

    p_elbow <- ElbowPlot(merged_sct, ndims = 30)
    ggsave("plots/ACP_elbow.png", p_elbow, width = 7, height = 5)

    # ── Harmony Integration ───────────────────────────────────────────────────
    integrated <- RunHarmony(merged_sct, group.by.vars = "sample",
                             reduction = "pca", reduction.save = "harmony",
                             verbose = FALSE)

    # ── UMAP + Clustering ─────────────────────────────────────────────────────
    integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:20,
                          reduction.name = "umap.harmony")
    integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:20)
    integrated <- FindClusters(integrated, resolution = 0.5, algorithm = 4)

    # ── UMAP Plot ─────────────────────────────────────────────────────────────
    p_umap <- DimPlot(integrated, reduction = "umap.harmony",
                      group.by = c("sample", "seurat_clusters"),
                      ncol = 2)
    ggsave("plots/ACP_integrated_UMAP.png", p_umap, width = 14, height = 6)

    saveRDS(integrated, "ACP_Integrated.rds")
    message("ACP integration complete: ", ncol(integrated), " cells")
    REOF
    """
}
