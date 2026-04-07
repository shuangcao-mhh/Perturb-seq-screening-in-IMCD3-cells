# 1. load the library 
library(scater)
library(Seurat)
library(tidyverse)
library(dplyr)
library(SingleCellExperiment)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(dittoSeq)


# 2. Read cellranger output and create Seurat object
CTL1.data <- Read10X(data.dir = "/hpc/project/ag-schmidt-ott/Shuang/3.Materials_for_paper_publish/1.Felix_WT_mus_snseq/3.cellranger_output/raw_feature_bc_matrix_Ctrl1")
CTL2.data <- Read10X(data.dir = "/hpc/project/ag-schmidt-ott/Shuang/3.Materials_for_paper_publish/1.Felix_WT_mus_snseq/3.cellranger_output/raw_feature_bc_matrix_Ctrl2")
TOL1.data <- Read10X(data.dir = "/hpc/project/ag-schmidt-ott/Shuang/3.Materials_for_paper_publish/1.Felix_WT_mus_snseq/3.cellranger_output/raw_feature_bc_matrix_Tol1")
TOL2.data <- Read10X(data.dir = "/hpc/project/ag-schmidt-ott/Shuang/3.Materials_for_paper_publish/1.Felix_WT_mus_snseq/3.cellranger_output/raw_feature_bc_matrix_Tol2")

CTL1 <- CreateSeuratObject(counts = CTL1.data, project = "CTL1", min.cells = 3, min.features = 500)
CTL2 <- CreateSeuratObject(counts = CTL2.data, project = "CTL2", min.cells = 3, min.features = 500)
TOL1 <- CreateSeuratObject(counts = TOL1.data, project = "TOL1", min.cells = 3, min.features = 500)
TOL2 <- CreateSeuratObject(counts = TOL2.data, project = "TOL2", min.cells = 3, min.features = 500)

CTL1[["percent.mt"]] <- PercentageFeatureSet(CTL1, pattern = "^mt-")
CTL2[["percent.mt"]] <- PercentageFeatureSet(CTL2, pattern = "^mt-")
TOL1[["percent.mt"]] <- PercentageFeatureSet(TOL1, pattern = "^mt-")
TOL2[["percent.mt"]] <- PercentageFeatureSet(TOL2, pattern = "^mt-")

p1 <- VlnPlot(CTL1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p2 <- VlnPlot(CTL2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p3 <- VlnPlot(TOL1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
p4 <- VlnPlot(TOL2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

p1|p2
p3|p4


# 3. Primarily filtration 
CT1 <- subset(CTL1, subset = nFeature_RNA > 600  & nFeature_RNA < 8000 & percent.mt < 2.5)
CT2 <- subset(CTL2, subset = nFeature_RNA > 600  & nFeature_RNA < 8000 & percent.mt < 2.5)
TO1 <- subset(TOL1, subset = nFeature_RNA > 600  & nFeature_RNA < 8000 & percent.mt < 2.5)
TO2 <- subset(TOL2, subset = nFeature_RNA > 600  & nFeature_RNA < 8000 & percent.mt < 2.5)

# 4. Merge data to one object 
Merge <- merge(CT1, c(CT2, TO1, TO2), add.cell.ids=c("CT1","CT2","TO1","TO2"))
table(Idents(Merge))

# 5. Integration 
Merge.list <- SplitObject(Merge, split.by = "orig.ident")

Merge.list <- lapply(X = Merge.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

Merge.anchors <- FindIntegrationAnchors(object.list = Merge.list, dims = 1:30)
Merge.combined <- IntegrateData(anchorset = Merge.anchors, dims = 1:30)

# 6. Integrated data process 
DefaultAssay(Merge.combined) <- "integrated"
# Run the standard workflow for visualization and clustering
Merge.combined <- ScaleData(Merge.combined, verbose = FALSE)
Merge.combined <- RunPCA(Merge.combined, npcs = 30, verbose = FALSE)
# uMAP and Clustering
Merge.combined <- RunUMAP(Merge.combined, reduction = "pca", dims = 1:30, seed.use = 123)
Merge.combined <- FindNeighbors(Merge.combined, reduction = "pca", dims = 1:30)
Merge.combined <- FindClusters(Merge.combined, resolution = 0.5)
# Visualization
px <- DimPlot(Merge.combined, reduction = "umap", group.by = "orig.ident")
py<- DimPlot(Merge.combined, reduction = "umap", label = F)
plot_grid(py)

DimPlot(Merge.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE)


#Initial Clusters - Broad Cell Types----
Merge.combined <- JoinLayers(Merge.combined)

cluster.markers <- list()

for (i in 0:25) {
  cluster.markers[[paste0("cluster", i, ".markers")]] <- FindConservedMarkers(
    Merge.combined,
    ident.1 = i,
    grouping.var = "orig.ident",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = TRUE
  )
}

#Save all cell type marker to csv document
all_markers <- list()

for (i in 0:25) {
  obj_name <- paste0("cluster", i, ".markers")
  df <- get(obj_name)
  
  df$gene <- rownames(df)   
  rownames(df) <- NULL      
  
  df$cluster <- i
  all_markers[[i + 1]] <- df
}

all_markers <- do.call(rbind, all_markers)

write.csv(all_markers, "all_clusters_markers.csv", row.names = FALSE)



#Rename the clusters according to the markers
Idents(Merge.combined) <- "seurat_clusters"
Merge.combined <- RenameIdents(Merge.combined, 
                               `0` = "PT", 
                               `1` = "PT", 
                               `2` = "TAL", 
                               `3` = "PT", 
                               `4` = "PT", 
                               `5` = "PT", 
                               `6` = "Endothelium", 
                               `7` = "DCT", 
                               `8` = "TAL", 
                               `9` = "CD-PC", 
                               `10` = "CNT", 
                               `11` = "TAL", 
                               `12` = "TL", 
                               `13` = "TAL/DCT", 
                               `14` = "PEC", 
                               `15` = "CD-IC",
                               `16` = "CD-PC",
                               `17` = "Interst.Cells", 
                               `18` = "CD-PC", 
                               `19` = "Podocytes", 
                               `20` = "TAL", 
                               `21` = "DME", 
                               `22` = "CD-PC",
                               `23` = "Immune.Cells", 
                               `24` = "Endothelium",
                               `25` = "TAL")

DimPlot(Merge.combined, reduction = "umap", pt.size = 0, label = TRUE)
DimPlot(Merge.combined, reduction = "umap", pt.size = 0, label = TRUE, group.by = "seurat_clusters")
Merge.combined$cell_type <- Idents(Merge.combined)
unique(Merge.combined@meta.data$cell_type)

#Conserved markers for Broad Cell Types----
DefaultAssay(Merge.combined) <- "RNA"
my_levels <- c(
               "Podocytes",
               "PEC",
               "PT",              
               "TL",             
               "TAL",              
               "TAL/DCT",             
               "DCT",
               "CNT",
               "CD-PC",         
               "CD-IC",  
               "DME",
               "Endothelium",     
               "Interst.Cells",  
               "Immune.Cells"    
)

Merge.combined$cell_type <- factor(Merge.combined$cell_type, levels = my_levels)
Idents(Merge.combined) <- "cell_type"

DotPlot(Merge.combined, features = c("Wt1", "Nphs2",                      # Podocytes
                               "Cldn1","Ncam1",                           # PEC
                               "Sord", "Cyp2e1", "Lrp2",                  # PT
                               "Bst1", "Slc4a11","Sptssb",                #Thin Limb
                               "Slc12a1","Umod",                          # TAL
                               "Slc12a3", "Calb1",                        # DCT & CNT
                               "Aqp2", "Aqp3", "Fxyd4",                   # CCD
                               "Aqp4", "Bcat1", "Tspan6", "Plat","Nupr1", # MCD
                               "Upk3a",                                   # DME
                               "Atp6v1g3", "Atp6v0d2",                    # CD-IC
                               "Pecam1", "Flt1","Cdh5",                   # Endothelium
                               "Pdgfrb", "Dcn","Fn1",                     # Interst. Cells
                               "H2-Ab1", "H2-Aa"), 
                               cols = c("white", "blue")) + RotatedAxis()


# Subset CD-PC  
levels(Idents(Merge.combined))
table(Idents(Merge.combined))

Idents(Merge.combined) <- "cell_type"
CDPC <- subset(Merge.combined, idents = 'CD-PC')
CDPC <- FindNeighbors(CDPC, dims = 1:30,reduction = "pca",k.param = 20)
CDPC <- FindClusters(CDPC, graph.name = "integrated_snn")

DimPlot(CDPC, reduction = "umap", label = T, seed = (123))

# Check the CD-PC sub-cluster markers 
CDPC <- JoinLayers(CDPC)

markers_cdpc <- FindAllMarkers(
  CDPC,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
write.csv(markers_cdpc, "all_cdpc_clusters_markers.csv", row.names = FALSE)


sessionInfo()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
> sessionInfo()
R version 4.4.2 (2024-10-31)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Etc/UTC
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] metap_1.13                  multtest_2.62.0             dittoSeq_1.18.0             cowplot_1.2.0              
 [5] future_1.70.0               RColorBrewer_1.1-3          DESeq2_1.46.0               pheatmap_1.0.13            
 [9] lubridate_1.9.5             forcats_1.0.1               stringr_1.6.0               dplyr_1.2.0                
[13] purrr_1.2.1                 readr_2.2.0                 tidyr_1.3.2                 tibble_3.3.1               
[17] tidyverse_2.0.0             Seurat_5.4.0                SeuratObject_5.3.0          sp_2.2-1                   
[21] scater_1.34.1               ggplot2_4.0.2               scuttle_1.16.0              SingleCellExperiment_1.28.1
[25] SummarizedExperiment_1.36.0 Biobase_2.66.0              GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
[29] IRanges_2.40.1              S4Vectors_0.44.0            BiocGenerics_0.52.0         MatrixGenerics_1.18.1      
[33] matrixStats_1.5.0          

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.23        splines_4.4.2           later_1.4.8             R.oo_1.27.1            
  [5] polyclip_1.10-7         fastDummies_1.7.5       lifecycle_1.0.5         Rdpack_2.6.6           
  [9] globals_0.19.1          lattice_0.22-6          MASS_7.3-65             magrittr_2.0.4         
 [13] limma_3.62.2            plotly_4.12.0           plotrix_3.8-14          qqconf_1.3.2           
 [17] httpuv_1.6.17           otel_0.2.0              sn_2.1.3                sctransform_0.4.3      
 [21] spam_2.11-3             sessioninfo_1.2.3       spatstat.sparse_3.1-0   reticulate_1.45.0      
 [25] pbapply_1.7-4           multcomp_1.4-30         abind_1.4-8             zlibbioc_1.52.0        
 [29] Rtsne_0.17              presto_1.0.0            R.utils_2.13.0          TH.data_1.1-5          
 [33] sandwich_3.1-1          GenomeInfoDbData_1.2.13 ggrepel_0.9.8           irlba_2.3.7            
 [37] listenv_0.10.1          spatstat.utils_3.2-2    TFisher_0.2.0           goftest_1.2-3          
 [41] RSpectra_0.16-2         spatstat.random_3.4-5   fitdistrplus_1.2-6      parallelly_1.46.1      
 [45] codetools_0.2-20        DelayedArray_0.32.0     tidyselect_1.2.1        UCSC.utils_1.2.0       
 [49] farver_2.1.2            ScaledMatrix_1.14.0     viridis_0.6.5           spatstat.explore_3.8-0 
 [53] mathjaxr_2.0-0          jsonlite_2.0.0          BiocNeighbors_2.0.1     progressr_0.19.0       
 [57] ggridges_0.5.7          survival_3.8-6          tools_4.4.2             ica_1.0-3              
 [61] Rcpp_1.1.1              glue_1.8.0              mnormt_2.1.2            gridExtra_2.3          
 [65] SparseArray_1.6.2       numDeriv_2016.8-1.1     withr_3.0.2             BiocManager_1.30.27    
 [69] fastmap_1.2.0           digest_0.6.39           rsvd_1.0.5              timechange_0.4.0       
 [73] R6_2.6.1                mime_0.13               colorspace_2.1-2        scattermore_1.2        
 [77] tensor_1.5.1            spatstat.data_3.1-9     R.methodsS3_1.8.2       generics_0.1.4         
 [81] data.table_1.18.2.1     httr_1.4.8              htmlwidgets_1.6.4       S4Arrays_1.6.0         
 [85] uwot_0.2.4              pkgconfig_2.0.3         gtable_0.3.6            lmtest_0.9-40          
 [89] S7_0.2.1                XVector_0.46.0          htmltools_0.5.9         dotCall64_1.2          
 [93] scales_1.4.0            png_0.1-9               spatstat.univar_3.1-7   rstudioapi_0.18.0      
 [97] tzdb_0.5.0              reshape2_1.4.5          nlme_3.1-169            zoo_1.8-15             
[101] KernSmooth_2.23-26      parallel_4.4.2          miniUI_0.1.2            vipor_0.4.7            
[105] ggrastr_1.0.2           pillar_1.11.1           grid_4.4.2              vctrs_0.7.2            
[109] RANN_2.6.2              promises_1.5.0          BiocSingular_1.22.0     beachmat_2.22.0        
[113] xtable_1.8-8            cluster_2.1.8.2         beeswarm_0.4.0          mvtnorm_1.3-6          
[117] cli_3.6.5               locfit_1.5-9.12         compiler_4.4.2          rlang_1.1.7            
[121] crayon_1.5.3            mutoss_0.1-14           future.apply_1.20.2     labeling_0.4.3         
[125] plyr_1.8.9              ggbeeswarm_0.7.3        stringi_1.8.7           viridisLite_0.4.3      
[129] deldir_2.0-4            BiocParallel_1.40.2     lazyeval_0.2.2          spatstat.geom_3.7-3    
[133] Matrix_1.7-5            RcppHNSW_0.6.0          hms_1.1.4               patchwork_1.3.2        
[137] statmod_1.5.1           shiny_1.13.0            rbibutils_2.4.1         ROCR_1.0-12            
[141] igraph_2.2.2








