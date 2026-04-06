# Analysis of snRNA-seq Data for [Interleukin Enhancer Binding Factor 2 (Ilf2) and Kidney Epithelial Stress Resilience]

[![BioRxiv](https://img.shields.io/badge/bioRxiv-10.1101/2025.08.12.667599-red)](https://www.biorxiv.org/content/10.1101/2025.08.12.667599v3)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete computational pipeline for the bioinformatic analysis presented in the manuscript. We have organized the scripts sequentially to ensure full reproducibility of our findings.

---

## 📂 Repository Structure & Workflow

To reproduce the analysis, please follow the scripts in the numerical order prefix:

| Step | Script Name | Description | Key Output |
| :--- | :--- | :--- | :--- |
| **01** | `01_mouse_kidney_snSeq_Preprocessing.R` | Data loading, MT-filtering. | Cleaned Seurat Object |
| **02** | `02_Integration_Normalization.R` | Batch effect correction and SCTransform. | Integrated Data |
| **03** | `03_Clustering_UMAP.R` | **Generation of Figure 1a.** UMAP & Clustering. | `Fig1a_UMAP.pdf` |
| **04** | `04_Cell_Annotation.R` | Marker gene identification and cluster naming. | Annotated Clusters |
| **05** | `05_Differential_Expression.R` | DEG analysis across conditions/groups. | `Table_S1_DEGs.csv` |

---

## 📍 Key Analysis: Figure 1a (snRNA-seq UMAP)

Reviewers and readers interested specifically in **Figure 1a** can find the implementation in:
👉 `scripts/03_Clustering_UMAP.R`

**Key Parameters used:**
* **Normalization:** SCTransform
* **PCs:** 1-30
* **Resolution:** 0.8
* **Markers Used for Annotation:** *Marker1, Marker2, Marker3...*

---

## 📈 Quality Control (QC) Summary

We provide detailed QC reports for all analyzed datasets in the `/results/QC_plots/` directory. These include:
* **Violin Plots:** nFeature, nCount, and percent.mt per sample.
* **Scatter Plots:** Correlation between nCount and nFeature.
* **Post-filtering stats:** Summary table of cell counts before and after QC.

---

## 💻 Setup & Requirements

### Dependencies
* **R version:** 4.2.0+
* **Core Packages:** `Seurat (v4.x)`, `tidyverse`, `DoubletFinder`, `harmony`

### How to use
1. Clone the repo: `git clone https://github.com/yourusername/reponame.git`
2. Place raw matrix files in `./data/raw/`
3. Run scripts in order: `Rscript scripts/01_QC_Preprocessing.R`

---

## 📧 Contact
For questions regarding the code or data, please open an issue or contact: **Shuang Cao** (cao.shuang@mh-hannover.de)
