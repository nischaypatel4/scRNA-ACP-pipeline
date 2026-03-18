/*
================================================================================
  Module: ACP_PERIO_INTEGRATION
  Purpose : Integrate ACP and Periodontium datasets using three methods:
            Harmony, CCA-LogNormalize, and CCA-SCTransform
  Input   : ACP and Periodontium merged + integrated Seurat objects
  Output  : Three integrated Seurat objects (one per method)
================================================================================
*/

process ACP_PERIO_INTEGRATION {
    publishDir "${params.outdir}/ACP_PERIO_Integrated", mode: 'copy'

    input:
    path acp_merged_raw
    path acp_integrated
    path perio_merged_raw
    path perio_integrated

    output:
    path "ACP_PERIO_Integrated_Harmony.rds",    emit: harmony
    path "ACP_PERIO_Integrated_CCA_Log.rds",    emit: cca_log
    path "ACP_PERIO_Integrated_CCA_SCT.rds",    emit: cca_sct
    path "ACP_PERIO_Merged_SCT.rds"
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

    acp_raw   <- readRDS("${acp_merged_raw}")
    acp_int   <- readRDS("${acp_integrated}")
    perio_raw <- readRDS("${perio_merged_raw}")
    perio_int <- readRDS("${perio_integrated}")

    # ── Method 1: Harmony Integration ────────────────────────────────────────
    merged_raw <- merge(acp_raw, perio_raw,
                        add.cell.ids = c("ACP", "PERIO"))
    merged_sct <- SCTransform(merged_raw, vars.to.regress = "percent.mt",
                              verbose = FALSE)
    saveRDS(merged_sct, "ACP_PERIO_Merged_SCT.rds")

    merged_sct <- RunPCA(merged_sct, npcs = 30, verbose = FALSE)
    harmony_int <- RunHarmony(merged_sct, group.by.vars = "orig.ident",
                              reduction = "pca", reduction.save = "harmony",
                              verbose = FALSE)
    harmony_int <- RunUMAP(harmony_int, reduction = "harmony", dims = 1:20,
                           reduction.name = "umap.harmony")
    harmony_int <- FindNeighbors(harmony_int, reduction = "harmony", dims = 1:20)
    harmony_int <- FindClusters(harmony_int, resolution = 0.5, algorithm = 4)

    p1 <- DimPlot(harmony_int, reduction = "umap.harmony",
                  group.by = "orig.ident") + ggtitle("Harmony: by Dataset")
    ggsave("plots/ACP_PERIO_Harmony_UMAP.png", p1, width = 8, height = 6)
    saveRDS(harmony_int, "ACP_PERIO_Integrated_Harmony.rds")

    # ── Method 2: CCA LogNormalize ────────────────────────────────────────────
    acp_ln  <- NormalizeData(acp_raw)  |> FindVariableFeatures()
    perio_ln <- NormalizeData(perio_raw) |> FindVariableFeatures()

    anchors_log <- FindIntegrationAnchors(
        object.list = list(acp_ln, perio_ln),
        dims = 1:20, reduction = "cca"
    )
    cca_log <- IntegrateData(anchorset = anchors_log, dims = 1:20)
    cca_log <- ScaleData(cca_log) |>
               RunPCA(npcs = 30, verbose = FALSE) |>
               RunUMAP(dims = 1:20) |>
               FindNeighbors(dims = 1:20) |>
               FindClusters(resolution = 0.5)

    p2 <- DimPlot(cca_log, group.by = "orig.ident") +
          ggtitle("CCA-LogNorm: by Dataset")
    ggsave("plots/ACP_PERIO_CCA_Log_UMAP.png", p2, width = 8, height = 6)
    saveRDS(cca_log, "ACP_PERIO_Integrated_CCA_Log.rds")

    # ── Method 3: CCA SCTransform ─────────────────────────────────────────────
    acp_sct  <- SCTransform(acp_raw,   vars.to.regress = "percent.mt", verbose = FALSE)
    perio_sct <- SCTransform(perio_raw, vars.to.regress = "percent.mt", verbose = FALSE)

    anchors_sct <- FindIntegrationAnchors(
        object.list = list(acp_sct, perio_sct),
        normalization.method = "SCT", dims = 1:20
    )
    cca_sct <- IntegrateData(anchorset = anchors_sct,
                             normalization.method = "SCT", dims = 1:20)
    cca_sct <- RunPCA(cca_sct, npcs = 30, verbose = FALSE) |>
               RunUMAP(dims = 1:20) |>
               FindNeighbors(dims = 1:20) |>
               FindClusters(resolution = 0.5)

    p3 <- DimPlot(cca_sct, group.by = "orig.ident") +
          ggtitle("CCA-SCT: by Dataset")
    ggsave("plots/ACP_PERIO_CCA_SCT_UMAP.png", p3, width = 8, height = 6)
    saveRDS(cca_sct, "ACP_PERIO_Integrated_CCA_SCT.rds")

    message("ACP-Periodontium integration complete")
    REOF
    """
}
