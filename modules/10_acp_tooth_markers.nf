process ACP_TOOTH_MARKERS {
    publishDir "${params.outdir}/ACP_Tooth_Integration/markers", mode: 'copy'

    input:
    path acp_integrated
    path tooth_rds

    output:
    path "*.png"
    path "*.csv"

    script:
    """
    Rscript - <<'REOF'
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(dplyr)

    set.seed(${params.seed})

    acp_int <- readRDS("${acp_integrated}")
    tooth   <- readRDS("${tooth_rds}")

    # ── ACP UMAP by subcelltype ───────────────────────────────────────────────
    p_acp <- DimPlot(acp_int, reduction = "umap.harmony",
                     group.by = "subcelltype", label = TRUE) +
             ggtitle("ACP: Subcell types")
    ggsave("ACP_subcelltype_UMAP.png", p_acp, width = 9, height = 7)

    # ── Tooth UMAP by cell type ───────────────────────────────────────────────
    if ("celltype" %in% colnames(tooth@meta.data)) {
        p_tooth <- DimPlot(tooth, group.by = "celltype", label = TRUE) +
                   ggtitle("Tooth: Cell types")
        ggsave("Tooth_celltype_UMAP.png", p_tooth, width = 9, height = 7)
    }

    # ── WNT signaling markers (key pathway in ACP) ────────────────────────────
    wnt_markers <- c("CTNNB1", "LEF1", "AXIN2", "WNT5A", "FZD7", "MYC")
    p_wnt <- FeaturePlot(acp_int, features = wnt_markers,
                         reduction = "umap.harmony", ncol = 3)
    ggsave("ACP_WNT_markers.png", p_wnt, width = 14, height = 9)

    # ── Enamel knot / epithelial markers (shared ACP-tooth features) ──────────
    enamel_markers <- c("PITX2", "BMP4", "SHH", "FGF4", "AMBN", "AMELX")
    p_enamel <- FeaturePlot(acp_int, features = enamel_markers,
                            reduction = "umap.harmony", ncol = 3)
    ggsave("ACP_enamel_markers.png", p_enamel, width = 14, height = 9)

    # ── Dot plot: top markers across subcell types ────────────────────────────
    Idents(acp_int) <- "subcelltype"
    top_markers <- FindAllMarkers(acp_int, only.pos = TRUE,
                                  min.pct = 0.25, logfc.threshold = 0.3)
    write.csv(top_markers, "ACP_subcelltype_markers.csv", row.names = FALSE)

    top5 <- top_markers |> group_by(cluster) |> slice_max(avg_log2FC, n = 5)
    p_dot <- DotPlot(acp_int, features = unique(top5\$gene)) +
             RotatedAxis() + ggtitle("Top 5 markers per subcell type")
    ggsave("ACP_marker_dotplot.png", p_dot, width = 14, height = 6)

    message("Marker visualisation complete")
    REOF
    """
}