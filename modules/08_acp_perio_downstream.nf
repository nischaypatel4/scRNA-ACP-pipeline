process ACP_PERIO_DOWNSTREAM {
    publishDir "${params.outdir}/ACP_PERIO_Integrated/downstream", mode: 'copy'

    input:
    path harmony_rds

    output:
    path "*.csv"
    path "*.png"

    script:
    """
    Rscript - <<'REOF'
    library(Seurat)
    library(clusterProfiler)
    library(ReactomePA)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
    library(dplyr)

    set.seed(${params.seed})

    obj <- readRDS("${harmony_rds}")

    # ── Marker identification per cluster ─────────────────────────────────────
    Idents(obj) <- "seurat_clusters"
    markers <- FindAllMarkers(obj, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25)
    write.csv(markers, "cluster_markers.csv", row.names = FALSE)

    # ── GO Enrichment on top cluster markers ──────────────────────────────────
    top_markers <- markers |>
        group_by(cluster) |>
        slice_max(avg_log2FC, n = 100) |>
        pull(gene) |>
        unique()

    entrez_ids <- bitr(top_markers, fromType = "SYMBOL",
                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)

    go_bp <- enrichGO(gene          = entrez_ids\$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05)

    write.csv(as.data.frame(go_bp), "GO_BP_enrichment.csv", row.names = FALSE)

    p_go <- dotplot(go_bp, showCategory = 20) + ggtitle("GO Biological Process")
    ggsave("GO_BP_dotplot.png", p_go, width = 10, height = 10)

    # ── Reactome Pathway Enrichment ───────────────────────────────────────────
    reactome <- enrichPathway(gene         = entrez_ids\$ENTREZID,
                              pvalueCutoff = 0.05,
                              readable     = TRUE)
    write.csv(as.data.frame(reactome), "Reactome_enrichment.csv",
              row.names = FALSE)

    p_reactome <- dotplot(reactome, showCategory = 20) +
                  ggtitle("Reactome Pathways")
    ggsave("Reactome_dotplot.png", p_reactome, width = 10, height = 10)

    message("Downstream analysis complete")
    REOF
    """
}


/*
================================================================================
  Module: ACP_TOOTH_INTEGRATION
  Purpose : Integrate ACP and Tooth developmental datasets using Harmony and RPCA
  Input   : ACP merged raw, ACP integrated, Tooth processed objects
  Output  : Harmony and RPCA integrated Seurat objects
================================================================================
*/