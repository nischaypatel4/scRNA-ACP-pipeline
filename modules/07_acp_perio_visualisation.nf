process ACP_PERIO_VIZ {
    publishDir "${params.outdir}/ACP_PERIO_Integrated/plots", mode: 'copy'

    input:
    path harmony_rds
    path cca_sct_rds

    output:
    path "*.png"

    script:
    """
    Rscript - <<'REOF'
    library(Seurat)
    library(ggplot2)
    library(patchwork)
    library(RColorBrewer)

    set.seed(${params.seed})

    harmony <- readRDS("${harmony_rds}")
    cca_sct <- readRDS("${cca_sct_rds}")

    # ── UMAP by cluster, dataset, and cell type ───────────────────────────────
    p1 <- DimPlot(harmony, reduction = "umap.harmony",
                  group.by = "seurat_clusters", label = TRUE) +
          ggtitle("Harmony: Clusters")

    p2 <- DimPlot(harmony, reduction = "umap.harmony",
                  group.by = "orig.ident") +
          ggtitle("Harmony: Dataset")

    p3 <- DimPlot(cca_sct, group.by = "seurat_clusters", label = TRUE) +
          ggtitle("CCA-SCT: Clusters")

    p4 <- DimPlot(cca_sct, group.by = "orig.ident") +
          ggtitle("CCA-SCT: Dataset")

    combined <- (p1 | p2) / (p3 | p4)
    ggsave("ACP_PERIO_comparison_UMAP.png", combined, width = 16, height = 12)

    # ── Feature plots for key markers ─────────────────────────────────────────
    key_markers <- c("EPCAM", "VIM", "CDH1", "KRT5", "SOX2", "CTNNB1")
    p_feat <- FeaturePlot(harmony, features = key_markers,
                          reduction = "umap.harmony", ncol = 3)
    ggsave("ACP_PERIO_key_markers.png", p_feat, width = 15, height = 10)

    message("Visualisation complete")
    REOF
    """
}


/*
================================================================================
  Module: ACP_PERIO_DOWNSTREAM
  Purpose : Differential expression and GO/Reactome enrichment analysis
            between ACP and Periodontium epithelial populations
================================================================================
*/