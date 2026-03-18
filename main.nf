#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
================================================================================
  scRNA-seq Pipeline: Comparative Single-Cell Analysis of ACP and
  Tooth Developmental Programs
  Author : Nischay Patel | IIT Kanpur | Amgen Scholars Program, Tsinghua 2025
================================================================================
*/

include { ACP_QC            } from './modules/01_acp_qc'
include { ACP_INTEGRATION   } from './modules/02_acp_integration'
include { SCIMILARITY_PREP  } from './modules/03_scimilarity_prep'
include { PERIO_QC          } from './modules/04_perio_qc_integration'
include { TOOTH_QC          } from './modules/05_tooth_qc_integration'
include { ACP_PERIO_INTEGRATION  } from './modules/06_acp_perio_integration'
include { ACP_PERIO_VIZ     } from './modules/07_acp_perio_visualisation'
include { ACP_PERIO_DOWNSTREAM   } from './modules/08_acp_perio_downstream'
include { ACP_TOOTH_INTEGRATION  } from './modules/09_acp_tooth_integration'
include { ACP_TOOTH_MARKERS } from './modules/10_acp_tooth_markers'

// ── Parameters ───────────────────────────────────────────────────────────────
params.acp1_dir       = "${params.data_dir}/raw/acp1"
params.acp2_dir       = "${params.data_dir}/raw/acp2"
params.acp3_dir       = "${params.data_dir}/raw/acp3"
params.perio_dir      = "${params.data_dir}/ref_datasets/GSE161266"
params.tooth_dir      = "${params.data_dir}/ref_datasets/GSE184749"
params.scimilarity_model = "${params.data_dir}/models/model_v1.1"
params.data_dir       = "data"
params.outdir         = "results"

// ── Workflow ──────────────────────────────────────────────────────────────────
workflow {

    // ── Step 1: ACP Quality Control and Filtering ───────────────────────────
    acp_dirs_ch = Channel.of(
        [params.acp1_dir, "acp1"],
        [params.acp2_dir, "acp2"],
        [params.acp3_dir, "acp3"]
    )
    acp_filtered = ACP_QC(acp_dirs_ch)

    // ── Step 2: ACP Integration (Harmony + SCTransform) ─────────────────────
    acp_integrated = ACP_INTEGRATION(
        acp_filtered.filter { it[1] == "acp1" }.map { it[0] },
        acp_filtered.filter { it[1] == "acp2" }.map { it[0] },
        acp_filtered.filter { it[1] == "acp3" }.map { it[0] }
    )

    // ── Step 3: Scimilarity Input Preparation ───────────────────────────────
    SCIMILARITY_PREP(acp_integrated.acp_final_rds)

    // ── Step 4: Periodontium QC and Integration ──────────────────────────────
    perio_integrated = PERIO_QC(params.perio_dir)

    // ── Step 5: Tooth Dataset QC and Integration ─────────────────────────────
    tooth_integrated = TOOTH_QC(params.tooth_dir)

    // ── Step 6: ACP + Periodontium Integration ───────────────────────────────
    acp_perio_integrated = ACP_PERIO_INTEGRATION(
        acp_integrated.merged_raw,
        acp_integrated.integrated,
        perio_integrated.merged_raw,
        perio_integrated.integrated
    )

    // ── Step 7: ACP + Periodontium Visualisation ─────────────────────────────
    ACP_PERIO_VIZ(
        acp_perio_integrated.harmony,
        acp_perio_integrated.cca_sct
    )

    // ── Step 8: Downstream Analysis (DE + GO Enrichment) ────────────────────
    ACP_PERIO_DOWNSTREAM(acp_perio_integrated.harmony)

    // ── Step 9: ACP + Tooth Integration ──────────────────────────────────────
    acp_tooth_integrated = ACP_TOOTH_INTEGRATION(
        acp_integrated.merged_raw,
        acp_integrated.integrated,
        tooth_integrated.processed
    )

    // ── Step 10: Marker Expression in ACP-Tooth Object ───────────────────────
    ACP_TOOTH_MARKERS(
        acp_integrated.integrated,
        tooth_integrated.processed
    )
}
