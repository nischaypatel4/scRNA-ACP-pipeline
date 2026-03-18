# scRNA-seq Pipeline: Comparative Single-Cell Analysis of ACP and Tooth Developmental Programs

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.0-brightgreen)](https://www.nextflow.io/)
[![R](https://img.shields.io/badge/R-4.4.0-blue)](https://www.r-project.org/)

## Overview

This Nextflow DSL2 pipeline performs end-to-end single-cell RNA-seq analysis
comparing **Adamantinomatous Craniopharyngioma (ACP)** tumour cells with tooth
developmental reference datasets (Periodontium and Tooth) to identify
transcriptional links between ACP epithelial states and tooth developmental
cell populations.

The analysis was developed as part of the **Amgen Scholars Program at
Tsinghua University (2025)** under Prof. Yinqing Li.

---

## Pipeline Overview

```
Raw 10x Data (ACP1, ACP2, ACP3)
        │
        ▼
  01. ACP QC & Filtering
        │
        ▼
  02. ACP Integration (Harmony + SCTransform)
        │
        ├──► 03. Scimilarity Input Preparation
        │         (epithelial subset export → Python Scimilarity)
        │
Reference Datasets
        │
        ├──► 04. Periodontium QC & Integration (GSE161266)
        │
        └──► 05. Tooth QC & Integration (GSE184749)
                 │
                 ▼
  06. ACP + Periodontium Integration
      (Harmony | CCA-LogNorm | CCA-SCTransform)
        │
        ├──► 07. Visualisation (UMAP, feature plots)
        └──► 08. Downstream Analysis (DE + GO + Reactome)

  09. ACP + Tooth Integration (Harmony + RPCA)
        │
        └──► 10. Marker Expression & Transcriptional Links
```

---

## Key Findings

- Identified transcriptional links between ACP epithelial states (KE, WE, PE)
  and tooth enamel knot and reticulum populations
- Differential expression and GO enrichment revealed WNT signalling, cell
  adhesion, and morphogenesis genes as key regulators in ACP
- Scimilarity (128D triplet-loss embeddings trained on 23M+ cells) enabled
  robust cross-tissue cell type annotation

---

## Repository Structure

```
scRNA-ACP-pipeline/
├── main.nf                        # Top-level workflow entry point
├── nextflow.config                # Execution profiles and parameters
├── Dockerfile                     # Container definition
├── modules/
│   ├── 01_acp_qc.nf              # ACP quality control and filtering
│   ├── 02_acp_integration.nf     # ACP Harmony + SCTransform integration
│   ├── 03_scimilarity_prep.nf    # Epithelial subset export for Scimilarity
│   ├── 04_perio_qc_integration.nf  # Periodontium QC and integration
│   ├── 05_tooth_qc_integration.nf  # Tooth dataset QC and integration
│   ├── 06_acp_perio_integration.nf # ACP-Periodontium multi-method integration
│   ├── 07_acp_perio_visualisation.nf # UMAP and feature plots
│   ├── 08_acp_perio_downstream.nf  # DE, GO, and Reactome enrichment
│   ├── 09_acp_tooth_integration.nf # ACP-Tooth integration
│   └── 10_acp_tooth_markers.nf   # Marker expression visualisation
├── data/
│   ├── raw/                       # 10x CellRanger outputs (acp1, acp2, acp3)
│   ├── ref_datasets/              # GSE161266 (Perio), GSE184749 (Tooth)
│   ├── processed/                 # Intermediate Seurat objects (.rds)
│   └── scimilarity_input/         # CSV exports for Python Scimilarity
├── results/                       # Pipeline outputs (plots, tables)
└── docs/                          # Additional documentation
```

---

## Requirements

### Nextflow Pipeline (R-based steps)
- [Nextflow](https://www.nextflow.io/) ≥ 23.0
- Docker or Singularity (recommended)
- Or: R 4.4+ with packages listed in `Dockerfile`

### Scimilarity (Python, steps 03 + notebooks)
- Python 3.10
- conda environment: `environment.yaml`

```bash
conda env create -f environment.yaml
conda activate scimilarity
pip install scimilarity
```

---

## Quick Start

```bash
# Clone the repository
git clone https://github.com/nischaypatel4/scRNA-ACP-pipeline.git
cd scRNA-ACP-pipeline

# Place your data
# data/raw/acp1, data/raw/acp2, data/raw/acp3  (10x CellRanger output)
# data/ref_datasets/GSE161266/  (Periodontium)
# data/ref_datasets/GSE184749/  (Tooth)

# Run with Docker
nextflow run main.nf -profile docker --outdir results

# Run on HPC with Singularity
nextflow run main.nf -profile singularity,slurm --outdir results

# Resume after interruption
nextflow run main.nf -profile docker -resume
```

---

## Parameters

| Parameter    | Default    | Description                         |
|-------------|------------|-------------------------------------|
| `--data_dir` | `data`     | Root data directory                 |
| `--outdir`   | `results`  | Output directory                    |
| `--seed`     | `9867`     | Random seed for reproducibility     |

---

## Reference Datasets

| Dataset    | GEO Accession | Description                          |
|-----------|--------------|--------------------------------------|
| ACP       | sciadv_2023  | Adamantinomatous Craniopharyngioma (3 samples) |
| Perio     | GSE161266    | Human Periodontium scRNA-seq (5 samples)       |
| Tooth     | GSE184749    | Human Tooth Developmental Atlas                |

---

## Publications

- Single-cell RNA sequencing highlights intratumor heterogeneity in ACP.
  *Science Advances* (DOI: 10.1126/sciadv.adc8933)
- Defining human mesenchymal and epithelial heterogeneity.
  *eLife* (DOI: 10.7554/eLife.62810)
- A single-cell atlas of human teeth.
  *iScience* (DOI: 10.1016/j.isci.2021.102405)

---

## Author

**Nischay Patel**
B.Tech, Biological Sciences and Bioengineering, IIT Kanpur
Amgen Scholars Program, Tsinghua University, 2025

- [nischaypatel4@gmail.com](mailto:nischaypatel4@gmail.com)
- [nischaypatel.com](https://nischaypatel.com)
- [GitHub: nischaypatel4](https://github.com/nischaypatel4)
