process ACP_TOOTH_INTEGRATION {
    publishDir "${params.outdir}/ACP_Tooth_Integration", mode: 'copy'

    input:
    path acp_merged_raw
    path acp_integrated
    path tooth_rds

    output:
    path "ACP_Tooth_Integrated_Harmony_SCT.rds", emit: harmony
    path "ACP_Tooth_Integrated_RPCA_SCT.rds",    emit: rpca
    path "ACP_Tooth_merged.rds"
    path "plots/*.png"

    script:
    """
    mkdir -p plots

    Rscript - <<'REOF'
    library(Seurat)
    library(harmony)
    library(ggplot2)
    library(patchwork)

    set.seed(${params.seed})

    acp_raw <- readRDS("${acp_merged_raw}")
    tooth   <- readRDS("${tooth_rds}")

    acp_raw\$dataset  <- "ACP"
    tooth\$dataset    <- "Tooth"

    # ── Merge ─────────────────────────────────────────────────────────────────
    merged <- merge(acp_raw, tooth, add.cell.ids = c("ACP", "Tooth"))
    saveRDS(merged, "ACP_Tooth_merged.rds")

    # ── SCTransform ───────────────────────────────────────────────────────────
    merged_sct <- SCTransform(merged, vars.to.regress = "percent.mt",
                              verbose = FALSE)
    merged_sct <- RunPCA(merged_sct, npcs = 30, verbose = FALSE)

    # ── Harmony ───────────────────────────────────────────────────────────────
    harmony_int <- RunHarmony(merged_sct, group.by.vars = "dataset",
                              reduction = "pca", reduction.save = "harmony",
                              verbose = FALSE)
    harmony_int <- RunUMAP(harmony_int, reduction = "harmony", dims = 1:20,
                           reduction.name = "umap.harmony")
    harmony_int <- FindNeighbors(harmony_int, reduction = "harmony", dims = 1:20)
    harmony_int <- FindClusters(harmony_int, resolution = 0.5, algorithm = 4)

    p1 <- DimPlot(harmony_int, reduction = "umap.harmony",
                  group.by = "dataset") + ggtitle("Harmony: ACP vs Tooth")
    ggsave("plots/ACP_Tooth_Harmony_UMAP.png", p1, width = 8, height = 6)
    saveRDS(harmony_int, "ACP_Tooth_Integrated_Harmony_SCT.rds")

    # ── RPCA ──────────────────────────────────────────────────────────────────
    acp_sct   <- SCTransform(acp_raw, vars.to.regress = "percent.mt", verbose = FALSE)
    tooth_sct <- SCTransform(tooth,   vars.to.regress = "percent.mt", verbose = FALSE)
    acp_sct   <- RunPCA(acp_sct,   npcs = 30, verbose = FALSE)
    tooth_sct <- RunPCA(tooth_sct, npcs = 30, verbose = FALSE)

    anchors_rpca <- FindIntegrationAnchors(
        object.list = list(acp_sct, tooth_sct),
        normalization.method = "SCT",
        reduction = "rpca", dims = 1:20
    )
    rpca_int <- IntegrateData(anchorset = anchors_rpca,
                              normalization.method = "SCT", dims = 1:20)
    rpca_int <- RunPCA(rpca_int, npcs = 30, verbose = FALSE) |>
                RunUMAP(dims = 1:20) |>
                FindNeighbors(dims = 1:20) |>
                FindClusters(resolution = 0.5)

    p2 <- DimPlot(rpca_int, group.by = "dataset") + ggtitle("RPCA: ACP vs Tooth")
    ggsave("plots/ACP_Tooth_RPCA_UMAP.png", p2, width = 8, height = 6)
    saveRDS(rpca_int, "ACP_Tooth_Integrated_RPCA_SCT.rds")

    message("ACP-Tooth integration complete")
    REOF
    """
}


/*
================================================================================
  Module: ACP_TOOTH_MARKERS
  Purpose : Visualise key marker gene expression in ACP and Tooth integrated object
            to identify transcriptional links between ACP epithelial states
            and tooth developmental cell populations
================================================================================
*/