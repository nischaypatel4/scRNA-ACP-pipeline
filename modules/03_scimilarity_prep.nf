/*
================================================================================
  Module: SCIMILARITY_PREP
  Purpose : Extract epithelial cell subsets and export count matrices + metadata
            as input for the Scimilarity deep-learning cell type annotation tool
  Input   : Integrated ACP Seurat object
  Output  : Raw count CSVs and metadata CSVs per ACP sample per epithelial type
================================================================================
*/

process SCIMILARITY_PREP {
    publishDir "${params.outdir}/Scimilarity", mode: 'copy'

    input:
    path acp_integrated_rds

    output:
    path "*.csv"

    script:
    """
    Rscript - <<'REOF'
    library(Seurat)
    library(ggplot2)
    library(Matrix)

    set.seed(${params.seed})

    acp <- readRDS("${acp_integrated_rds}")

    # ── Normalise and cluster to identify epithelial populations ──────────────
    acp_swf <- NormalizeData(acp)
    acp_swf <- FindVariableFeatures(acp_swf)
    acp_swf <- ScaleData(acp_swf)
    acp_swf <- RunPCA(acp_swf, npcs = 30, verbose = FALSE)
    acp_swf <- RunUMAP(acp_swf, dims = 1:20)
    acp_swf <- FindNeighbors(acp_swf, dims = 1:20)
    acp_swf <- FindClusters(acp_swf, resolution = 0.3)

    # ── Epithelial subset labels expected from manual annotation ─────────────
    # (KE = Keratinised Epithelium, WE = Wet Epithelium,
    #  PE = Papillary Epithelium, CC = Ciliated Cells)
    epi_types <- c("Epi_KE", "Epi_WE", "Epi_PE", "Epi_CC")

    samples <- unique(acp\$sample)

    for (s in samples) {
        s_lower <- tolower(s)  # e.g. "ACP1" -> "acp1"

        for (epi in epi_types) {
            # Subset: cells from this sample AND this epithelial type
            cells <- WhichCells(acp_swf,
                                expression = sample == s & subcelltype == epi)

            if (length(cells) < 10) next  # Skip if too few cells

            sub <- subset(acp_swf, cells = cells)

            # Export raw counts
            raw_mat <- GetAssayData(sub, assay = "RNA", layer = "counts")
            write.csv(as.matrix(raw_mat),
                      paste0(s_lower, "_", tolower(epi), "_raw_counts.csv"))

            # Export metadata
            write.csv(sub@meta.data,
                      paste0(s_lower, "_", tolower(epi), "_metadata.csv"))

            message(s, " - ", epi, ": ", length(cells), " cells exported")
        }
    }
    REOF
    """
}
