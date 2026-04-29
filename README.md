# 📊 Analysis Pipeline for ILF2 Kidney Study

This repository contains the complete computational workflow for snRNA-seq, sc-CRISPR screening, and bulk RNA-seq analyses.
All scripts are organized sequentially and grouped by analysis type to ensure clarity and reproducibility.

---

# 🧠 Analysis Environment

This pipeline integrates **multiple programming environments**:

- **Bash**: alignment and preprocessing (e.g., STAR, Cellranger)
- **Python**: regulatory network inference and data processing (e.g., pySCENIC)
- **R**: downstream analysis, visualization, and statistics (e.g., Seurat, DESeq2)

---

# 🔁 Analysis Workflow

## 🧬 snRNA-seq Analysis

| Step | Script | Description | Output |
|------|--------|-------------|--------|
| 01 | `01_snRNA_processing` | QC, filtering, clustering, annotation (Seurat) (**Fig 1a**)| Annotated Seurat object |
| 02 | `02_snRNA_subcluster` | Sub-clustering (CD-PC, etc.) (**Fig5b**) | Sub-clusters |
| 03 | `03_map_IMCD3_to_reference` | Map IMCD3 cells to kidney reference (**Fig 1a**) | UMAP projection |
| 04 | `04_MCD_expression` | Gene expression in MCD cells (**Fig 1b, TF selection**) | Expression matrix |
| 05 | `05_SCENIC_analysis` | Regulon inference (SCENIC) (**Fig 1b, TF selection**) | TF networks |

---

## 🧪 sc-CRISPR Screen

| Step | Script | Description | Output |
|------|--------|----------------------|--------|
| 06 | `06_scCRISPR_QC` | QC and gRNA recovery (**Fig 1c, 2a**) | QC metrics |
| 07 | `07_scCRISPR_PCA_DEG` | PCA, Knockdown efficiency and DEGs (**Fig 1f, 2b**)| PCA_plot DEG tables |
| 08 | `08_scCRISPR_heatmap` | single-cell DEGs heatmap (**Fig 1f–g, 2b–c**) | sc-DEGs heatmap |
| 09 | `09_scCRISPR_ILF2_enrichment` | GO/KEGG enrichment (**Fig 2d**) | Pathways |
| 10 | `10_scCRISPR_ILF2_GSEA` | GSEA (**Fig 2e**) | Ranked pathways |

---

## 🧬 Bulk RNA-seq

| Step | Script | Description | Output |
|------|--------|-------------|--------|
| 11 | `11_bulk_IMCD3_expression` | Alignment (STAR), TPM (**Fig 1b, TF selection**) | Expression |
| 12 | `12_bulk_ILF2_DEG` | Alignment (STAR), TPM, DESeq2 (**Fig 3f–g**) | DEG results |
| 13 | `13_bulk_vs_sc_comparison` | Compare bulk vs sc (**Fig 3d–e**) | Correlation |
| 14 | `14_bulk_splicing_enrichment` | RNA splicing and enrichment(SpliceWiz) (**Fig 3h–k**) | Splicing results Pathways |
| 15 | `15_bulk_IRI_expression` | Alignment (STAR), TPM (**Fig 5g, 5j**) | Expression |
| 16 | `16_bulk_IRI_splicing` | Alternative splicing (SpliceWiz) (**Fig 5k, 5l**) | Splicing results |

---

# 🧰 Key Tools Used

- STAR
- Cellranger
- Seurat (snRNA-seq)
- Seurat / Mixscape (CRISPRi single-cell) 
- SCENIC (regulatory network)
- STAR (alignment)
- DESeq2 (DEG analysis)
- CombatSeq (batch correction)
- SpliceWiz (splicing)

---

# ⚠️ Notes
