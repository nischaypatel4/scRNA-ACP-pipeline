/*
================================================================================
  Module: ACP_QC
  Purpose : Quality control and filtering of ACP single-cell datasets
  Input   : 10x Genomics CellRanger output directory + sample name
  Output  : Filtered Seurat objects (.rds), QC plots
================================================================================
*/

process ACP_QC {
    tag        "${sample_id}"
    publishDir "${params.outdir}/ACP_QC_Filtering/${sample_id}", mode: 'copy'

    input:
    tuple path(data_dir), val(sample_id)

    output:
    tuple path("${sample_id}_filtered.rds"), val(sample_id), emit: filtered
    tuple path("${sample_id}_final.rds"),    val(sample_id), emit: final_rds
    path "plots/*.png"

    script:
    """
    mkdir -p plots

    Rscript - <<'REOF'
    library(Seurat)
    library(DropletUtils)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
    library(Matrix)
    library(scales)

    set.seed(${params.seed})

    # ── Load 10x data ─────────────────────────────────────────────────────────
    counts <- Read10X(data.dir = "${data_dir}")
    seurat_obj <- CreateSeuratObject(
        counts   = counts,
        project  = "${sample_id}",
        min.cells = 3,
        min.features = 200
    )

    # ── Calculate QC metrics ──────────────────────────────────────────────────
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

    # ── QC Violin Plot ────────────────────────────────────────────────────────
    p_vln <- VlnPlot(
        seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0
    )
    ggsave("plots/${sample_id}_qc_violin.png", p_vln, width = 12, height = 5)

    # ── Scatter plots ─────────────────────────────────────────────────────────
    p_scatter <- FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA") +
                 FeatureScatter(seurat_obj, "nCount_RNA", "percent.mt")
    ggsave("plots/${sample_id}_qc_scatter.png", p_scatter, width = 10, height = 5)

    # ── Barcode rank plot (knee plot) ─────────────────────────────────────────
    br_out <- barcodeRanks(counts)
    png("plots/${sample_id}_barcode_rank.png", width = 800, height = 600)
    plot(br_out\$rank, br_out\$total, log = "xy",
         xlab = "Barcode Rank", ylab = "Total UMI", main = "${sample_id}")
    abline(h = metadata(br_out)\$knee, col = "blue", lty = 2)
    abline(h = metadata(br_out)\$inflection, col = "red", lty = 2)
    dev.off()

    # ── Filter cells ─────────────────────────────────────────────────────────
    filtered <- subset(
        seurat_obj,
        subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 25
    )
    saveRDS(filtered, "${sample_id}_filtered.rds")

    # ── Doublet detection + final object ─────────────────────────────────────
    # Simplified: apply additional nCount filter as proxy for doublet removal
    final <- subset(filtered, subset = nCount_RNA < quantile(filtered\$nCount_RNA, 0.975))
    saveRDS(final, "${sample_id}_final.rds")

    message("${sample_id}: ", ncol(seurat_obj), " -> ", ncol(final), " cells after QC")
    REOF
    """
}
