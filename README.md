# 📊 Analysis Pipeline for ILF2 Kidney Study

This repository contains the complete computational workflow for snRNA-seq, sc-CRISPR screening, and bulk RNA-seq analyses.

All scripts are organized sequentially and grouped by analysis type to ensure clarity and reproducibility.

---

# 📁 Repository Structure
scripts/
├── snRNA_seq/
├── sc_CRISPR/
├── bulk_RNAseq/
---

# 🔁 Analysis Workflow

Run scripts in numerical order.

## 🧬 snRNA-seq Analysis

| Step | Script | Description | Output |
|------|--------|------------|--------|
| 01 | `01_snRNA_processing.R` | QC, filtering, clustering, annotation (Seurat) | Annotated Seurat object |
| 02 | `02_snRNA_subcluster.R` | Sub-clustering (CD-PC, etc.) | Sub-clusters |
| 03 | `03_map_IMCD3_to_reference.R` | Map IMCD3 cells to kidney reference (**Fig 1a**) | UMAP projection |
| 04 | `04_SCENIC_analysis.R` | Regulon inference (SCENIC) (**Fig 1b**) | TF networks |
| 05 | `05_MCD_expression.R` | Gene expression in MCD cells | Expression matrix |

---

## 🧪 sc-CRISPR Screen

| Step | Script | Description | Output |
|------|--------|------------|--------|
| 06 | `06_scCRISPR_QC.R` | QC and gRNA recovery (**Fig 1c, 2a**) | QC metrics |
| 07 | `07_scCRISPR_DEG.R` | Knockdown efficiency and DEGs | DEG tables |
| 08 | `08_scCRISPR_PCA_heatmap.R` | PCA and heatmap (**Fig 1f–g, 2b–c**) | PCA, heatmap |
| 09 | `09_scCRISPR_ILF2_enrichment.R` | GO/KEGG enrichment (**Fig 2d**) | Pathways |
| 10 | `10_scCRISPR_ILF2_GSEA.R` | GSEA (**Fig 2e**) | Ranked pathways |

---

## 🧬 Bulk RNA-seq

| Step | Script | Description | Output |
|------|--------|------------|--------|
| 11 | `11_bulk_ILF2_DEG.R` | Alignment (STAR), TPM, DESeq2 (**Fig 3f–g**) | DEG results |
| 12 | `12_bulk_vs_sc_comparison.R` | Compare bulk vs sc (**Fig 3d–e**) | Correlation |
| 13 | `13_bulk_splicing.R` | Alternative splicing (SpliceWiz) (**Fig 3h–j**) | Splicing results |
| 14 | `14_bulk_splicing_enrichment.R` | Splicing enrichment (**Fig 3k**) | Pathways |
| 15 | `15_IRI_expression.R` | IRI gene expression (**Fig 5g, 5j**) | Expression |
| 16 | `16_IRI_splicing.R` | IRI splicing (**Fig 5k, 5l**) | Splicing |

---

# 🧰 Tools Used

- Seurat (snRNA-seq)
- SCENIC (regulatory network)
- STAR (alignment)
- DESeq2 (DEG analysis)
- CombatSeq (batch correction)
- SpliceWiz (splicing)

---

# ⚠️ Notes

- **Figure 1a UMAP** → `03_map_IMCD3_to_reference.R`
- QC and preprocessing are included in:
  - snRNA-seq → Step 01  
  - sc-CRISPR → Step 06  
  - bulk RNA-seq → Step 11  

- Scripts should be run sequentially.

---
