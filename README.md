# 📊 Analysis of snRNA-seq, sc-CRISPR and Bulk RNA-seq Data for ILF2 Study

This repository contains the complete computational pipeline for the bioinformatic analyses presented in the manuscript.  
All scripts are organized sequentially to ensure clarity, reproducibility, and traceability between scripts and figures.

---

# 📁 Repository Structure & Workflow

To reproduce the analysis, please follow the scripts in numerical order:

| Step | Script Name | Description | Key Output |
|------|------------|-------------|------------|
| 01 | `1.Mouse_kidney_snSeq_processing.R` | Import Cell Ranger output into R and process using Seurat, including QC filtering (mitochondrial filtering), clustering, and initial cell-type annotation | Filtered and annotated Seurat object, QC metrics, marker genes |
| 02 | `2.Sub_cluster_of_mouse_kidney_seurat_object.R` | Subset major cell clusters and extract control (WT) samples; perform sub-clustering of CD-PC and other relevant cell types | Sub-clustered Seurat object, marker genes for sub-clusters |
| 03 | `3.Map_IMCD3_cells_to_mouse_kidney.R` | Map IMCD3 scRNA-seq data (from CRISPRi non-targeting control cells) onto WT mouse kidney reference to generate projected UMAP (Figure 1a) | UMAP projection with IMCD3 cells mapped to kidney reference |
| 04 | `4.SCENIC_analysis_on_WT_kidney_snSeq.R` | Perform SCENIC analysis on WT kidney snRNA-seq data to identify active regulons and cell-type-specific transcription factor networks | Regulon activity matrix and TF network results |
| 05 | `5.WT_kidney_MCD_gene_expression.R` | Calculate average gene expression within MCD sub-clusters for transcription factor selection | Normalized average gene expression (MCD cells) |
| 06 | `6.WT_IMCD3_cells_bulk_RNAseq.R` | Process bulk RNA-seq data from WT IMCD3 cells and quantify gene expression for transcription factor comparison | TPM-normalized gene expression matrix |
| 07 | `7.sc_CRISPR_1_QC_recovery_quantification.R` | Perform quality control and quantify gRNA recovery and abundance in sc-CRISPR screen data | QC metrics, gRNA counts, cell coverage statistics |
| 08 | `8.sc_CRISPR_2_Knockdown_effectiveness_and_DEGs_quantification.R` | Perform pseudo-bulk differential expression analysis to evaluate knockdown efficiency and identify DEGs | DEG tables, normalized expression values, statistical significance |
| 09 | `9.sc_CRISPR_3_PCA_and_DEGs_heatmap.R` | Generate PCA plots per gRNA and visualize DEGs using heatmaps at single-cell level | PCA plots and DEG heatmaps |
| 10 | `10.sc_CRISPR_4_ILF2_KD_Enrichment_analysis.R` | Perform GO and KEGG enrichment analysis using DEGs from ILF2 knockdown (pseudo-bulk sc-CRISPR data) | Functional enrichment results (GO, KEGG) |
| 11 | `11.sc_CRISPR_5_ILF2_KD_GSEA_analysis.R` | Perform GSEA using ranked gene lists based on differential expression after ILF2 knockdown | Ranked gene list and GSEA results |
| 12 | `12.ILF2_KD_bulk_1_gene_expression.R` | Process bulk RNA-seq data (STAR alignment), perform TPM normalization, batch correction (CombatSeq), and DESeq2 differential expression analysis | Normalized expression matrix, DEG tables, statistical results |
| 13 | `13.ILF2_KD_bulk_2_correlation_with_sc_screen.R` | Compare DEGs between bulk RNA-seq and sc-CRISPR datasets to assess consistency and robustness | Venn diagrams and logFC correlation plots |
| 14 | `14.ILF2_KD_bulk_3_alternative_splicing.R` | Perform alternative splicing analysis using SpliceWiz on ILF2 knockdown bulk RNA-seq data | Exon coverage plots and differential splicing results |
| 15 | `15.ILF2_KD_bulk_4_splicing_enrichment_analysis.R` | Perform GO and KEGG enrichment analysis on differentially spliced genes | Functional enrichment tables |
| 16 | `16.IRI_bulk_1_gene_expression.R` | Process IRI bulk RNA-seq data and perform TPM normalization across samples | Normalized gene expression matrix |
| 17 | `17.IRI_bulk_2_alternative_splicing.R` | Perform alternative splicing analysis on IRI samples using SpliceWiz | Differential splicing results and exon coverage plots |

---

# 📌 Notes on Analysis Tools

- snRNA-seq analysis: **Seurat**  
- Regulatory network inference: **SCENIC**  
- Alignment: **STAR**  
- Differential expression: **DESeq2**  
- Batch correction: **CombatSeq**  
- Alternative splicing: **SpliceWiz**

---

# ⚠️ Key Notes for Reproducibility

- **QC and preprocessing steps are explicitly included:**
  - snRNA-seq preprocessing → Step 01  
  - sc-CRISPR QC → Step 07  
  - bulk RNA-seq preprocessing → Step 12  

- The scripts are designed to be executed **sequentially**, where outputs from earlier steps are used as inputs for downstream analyses.

---

## 📧 Contact
For questions regarding the code or data, please open an issue or contact: **Shuang Cao** (cao.shuang@mh-hannover.de)
