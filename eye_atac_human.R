# eye_atac_2
# GSE183684
# Lyu et al. Cell Rep. 2021 PMID: 34788628

library(Seurat)
library(Signac)
library(ArchR)

addArchRGenome("hg38")
genomeAnnotation <- createGenomeAnnotation(genome = "hg38")

##############################################################################

RNA_d53 <- ReadMtx(mtx = "GSE183684_RAW/GSM5567525_d53_matrix.mtx.gz",
                   cells = "GSE183684_RAW/GSM5567525_d53_barcodes.tsv.gz",
                   features = "GSE183684_RAW/GSM5567525_d53_features.tsv.gz")

RNA_d53 <- CreateSeuratObject(counts = RNA_d53,
                              project = "RNA_d53",
                              min.cells = 3,
                              min.features = 200)

# RNA_d53
# An object of class Seurat 
# 25350 features across 5189 samples within 1 assay 
# Active assay: RNA (25350 features, 0 variable features)
# 2 layers present: counts, data

RNA_d53[["percent.mt"]] <- PercentageFeatureSet(RNA_d53, pattern = "^MT-")

VlnPlot(RNA_d53, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

RNA_d53 <- subset(RNA_d53, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 4)

# RNA_d53
# An object of class Seurat 
# 25350 features across 2199 samples within 1 assay 
# Active assay: RNA (25350 features, 0 variable features)
# 2 layers present: counts, data

RNA_d53 <- NormalizeData(RNA_d53, normalization.method = "LogNormalize", scale.factor = 10000)

RNA_d53 <- FindVariableFeatures(RNA_d53, selection.method = "vst", nfeatures = 2000)

all.genes_d53 <- rownames(RNA_d53)
RNA_d53 <- ScaleData(RNA_d53, features = all.genes_d53)

RNA_d53 <- RunPCA(RNA_d53, features = VariableFeatures(object = RNA_d53))

ElbowPlot(RNA_d53)

RNA_d53_PCA <- FindNeighbors(RNA_d53, dims = 1:10)
RNA_d53_res_1.0 <- FindClusters(RNA_d53_PCA, resolution = 1.0)

RNA_d53_res_1.0 <- RunUMAP(RNA_d53_res_1.0, dims = 1:5)
DimPlot(RNA_d53_res_1.0, reduction = "umap")

# saveRDS(RNA_d53_res_1.0, file = "RNA_d53.rds")

##############################################################################

RNA_d74 <- ReadMtx(mtx = "GSE183684_RAW/GSM5567527_d74_matrix.mtx.gz",
                   cells = "GSE183684_RAW/GSM5567527_d74_barcodes.tsv.gz",
                   features = "GSE183684_RAW/GSM5567527_d74_features.tsv.gz")

RNA_d74 <- CreateSeuratObject(counts = RNA_d74,
                              project = "RNA_d74",
                              min.cells = 3,
                              min.features = 200)

# RNA_d74
# An object of class Seurat 
# 25600 features across 11059 samples within 1 assay 
# Active assay: RNA (25600 features, 0 variable features)
# 2 layers present: counts, data

RNA_d74[["percent.mt"]] <- PercentageFeatureSet(RNA_d74, pattern = "^MT-")

VlnPlot(RNA_d74, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

RNA_d74 <- subset(RNA_d74, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 5)

# RNA_d74
# An object of class Seurat 
# 25600 features across 9038 samples within 1 assay 
# Active assay: RNA (25600 features, 0 variable features)
# 2 layers present: counts, data

RNA_d74 <- NormalizeData(RNA_d74, normalization.method = "LogNormalize", scale.factor = 10000)

RNA_d74 <- FindVariableFeatures(RNA_d74, selection.method = "vst", nfeatures = 2000)

all.genes_d74 <- rownames(RNA_d74)
RNA_d74 <- ScaleData(RNA_d74, features = all.genes_d74)

RNA_d74 <- RunPCA(RNA_d74, features = VariableFeatures(object = RNA_d74))

ElbowPlot(RNA_d74)

RNA_d74_PCA <- FindNeighbors(RNA_d74, dims = 1:10)
RNA_d74_res_0.5 <- FindClusters(RNA_d74_PCA, resolution = 0.5)

RNA_d74_res_0.5 <- RunUMAP(RNA_d74_res_0.5, dims = 1:10)
DimPlot(RNA_d74_res_0.5, reduction = "umap")

# saveRDS(RNA_d74_res_0.5, file = "RNA_d74.rds")

##############################################################################

RNA_d78 <- ReadMtx(mtx = "GSE183684_RAW/GSM5567528_d78_matrix.mtx.gz",
                   cells = "GSE183684_RAW/GSM5567528_d78_barcodes.tsv.gz",
                   features = "GSE183684_RAW/GSM5567528_d78_features.tsv.gz")

RNA_d78 <- CreateSeuratObject(counts = RNA_d78,
                              project = "RNA_d78",
                              min.cells = 3,
                              min.features = 200)

# RNA_d78
# An object of class Seurat 
# 24914 features across 8681 samples within 1 assay 
# Active assay: RNA (24914 features, 0 variable features)
# 2 layers present: counts, data

RNA_d78[["percent.mt"]] <- PercentageFeatureSet(RNA_d78, pattern = "^MT-")

VlnPlot(RNA_d78, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

RNA_d78 <- subset(RNA_d78, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 8)

# RNA_d78
# An object of class Seurat 
# 24914 features across 8262 samples within 1 assay 
# Active assay: RNA (24914 features, 0 variable features)
# 2 layers present: counts, data

RNA_d78 <- NormalizeData(RNA_d78, normalization.method = "LogNormalize", scale.factor = 10000)

RNA_d78 <- FindVariableFeatures(RNA_d78, selection.method = "vst", nfeatures = 2000)

all.genes_d78 <- rownames(RNA_d78)
RNA_d78 <- ScaleData(RNA_d78, features = all.genes_d78)

RNA_d78 <- RunPCA(RNA_d78, features = VariableFeatures(object = RNA_d78))

ElbowPlot(RNA_d78)

RNA_d78_PCA <- FindNeighbors(RNA_d78, dims = 1:15)
RNA_d78_res_0.5 <- FindClusters(RNA_d78_PCA, resolution = 0.5)

RNA_d78_res_0.5 <- RunUMAP(RNA_d78_res_0.5, dims = 1:15)
DimPlot(RNA_d78_res_0.5, reduction = "umap")

# saveRDS(RNA_d78_res_0.5, file = "RNA_d78.rds")

##############################################################################

RNA_d113 <- ReadMtx(mtx = "GSE183684_RAW/GSM5567529_d113_matrix.mtx.gz",
                   cells = "GSE183684_RAW/GSM5567529_d113_barcodes.tsv.gz",
                   features = "GSE183684_RAW/GSM5567529_d113_features.tsv.gz")

RNA_d113 <- CreateSeuratObject(counts = RNA_d113,
                              project = "RNA_d113",
                              min.cells = 3,
                              min.features = 200)

# RNA_d113
# An object of class Seurat 
# 23315 features across 4365 samples within 1 assay 
# Active assay: RNA (23315 features, 0 variable features)
# 2 layers present: counts, data

RNA_d113[["percent.mt"]] <- PercentageFeatureSet(RNA_d113, pattern = "^MT-")

VlnPlot(RNA_d113, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

RNA_d113 <- subset(RNA_d113, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 4)

# RNA_d113
# An object of class Seurat 
# 23315 features across 3753 samples within 1 assay 
# Active assay: RNA (23315 features, 0 variable features)
# 2 layers present: counts, data

RNA_d113 <- NormalizeData(RNA_d113, normalization.method = "LogNormalize", scale.factor = 10000)

RNA_d113 <- FindVariableFeatures(RNA_d113, selection.method = "vst", nfeatures = 2000)

all.genes_d113 <- rownames(RNA_d113)
RNA_d113 <- ScaleData(RNA_d113, features = all.genes_d113)

RNA_d113 <- RunPCA(RNA_d113, features = VariableFeatures(object = RNA_d113))

ElbowPlot(RNA_d113)

RNA_d113_PCA <- FindNeighbors(RNA_d113, dims = 1:15)
RNA_d113_res_0.5 <- FindClusters(RNA_d113_PCA, resolution = 0.5)

RNA_d113_res_0.5 <- RunUMAP(RNA_d113_res_0.5, dims = 1:15)
DimPlot(RNA_d113_res_0.5, reduction = "umap")

# saveRDS(RNA_d113_res_0.5, file = "RNA_d113.rds")

##############################################################################

RNA_d132 <- ReadMtx(mtx = "GSE183684_RAW/GSM5567530_d132_matrix.mtx.gz",
                    cells = "GSE183684_RAW/GSM5567530_d132_barcodes.tsv.gz",
                    features = "GSE183684_RAW/GSM5567530_d132_features.tsv.gz")

RNA_d132 <- CreateSeuratObject(counts = RNA_d132,
                               project = "RNA_d132",
                               min.cells = 3,
                               min.features = 200)

# RNA_d132
# An object of class Seurat 
# 22953 features across 5018 samples within 1 assay 
# Active assay: RNA (22953 features, 0 variable features)
# 2 layers present: counts, data

RNA_d132[["percent.mt"]] <- PercentageFeatureSet(RNA_d132, pattern = "^MT-")

VlnPlot(RNA_d132, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

RNA_d132 <- subset(RNA_d132, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 4)

# RNA_d132
# An object of class Seurat 
# 22953 features across 4629 samples within 1 assay 
# Active assay: RNA (22953 features, 0 variable features)
# 2 layers present: counts, data

RNA_d132 <- NormalizeData(RNA_d132, normalization.method = "LogNormalize", scale.factor = 10000)

RNA_d132 <- FindVariableFeatures(RNA_d132, selection.method = "vst", nfeatures = 2000)

all.genes_d132 <- rownames(RNA_d132)
RNA_d132 <- ScaleData(RNA_d132, features = all.genes_d132)

RNA_d132 <- RunPCA(RNA_d132, features = VariableFeatures(object = RNA_d132))

ElbowPlot(RNA_d132)

RNA_d132_PCA <- FindNeighbors(RNA_d132, dims = 1:10)
RNA_d132_res_0.5 <- FindClusters(RNA_d132_PCA, resolution = 0.5)

RNA_d132_res_0.5 <- RunUMAP(RNA_d132_res_0.5, dims = 1:10)
DimPlot(RNA_d132_res_0.5, reduction = "umap")

# saveRDS(RNA_d132_res_0.5, file = "RNA_d132.rds")

##############################################################################

# inputFile_d53 <- c("d53"="GSE183684_RAW/GSM5567517_d53_fragments.tsv.gz")
# 
# addArchRThreads(1)
# 
# ArrowFile_d53 <- createArrowFiles(
#   inputFiles = inputFile_d53,
#   sampleNames = names(inputFile_d53),
#   minTSS = 0,
#   minFrags = 1000, 
#   addTileMat = TRUE,
#   addGeneScoreMat = TRUE
# )
# 
# doubScores <- addDoubletScores(
#   input = ArrowFile_d53,
#   k = 10, 
#   knnMethod = "UMAP", 
#   LSIMethod = 1
# )
# 
# proj_d53 <- ArchRProject(
#   ArrowFiles = ArrowFile_d53, 
#   outputDirectory = "d53",
#   copyArrows = TRUE
# )
# 
# proj_d53 <- filterDoublets(ArchRProj = proj_d53)
# 
# proj_d53 <- saveArchRProject(ArchRProj = proj_d53)

##############################################################################

# inputFile_d74 <- c("d74"="GSE183684_RAW/GSM5567519_d74_fragments.tsv.gz")
# 
# addArchRThreads(1)
# 
# ArrowFile_d74 <- createArrowFiles(
#   inputFiles = inputFile_d74,
#   sampleNames = names(inputFile_d74),
#   minTSS = 0,
#   minFrags = 1000, 
#   addTileMat = TRUE,
#   addGeneScoreMat = TRUE
# )
# 
# doubScores <- addDoubletScores(
#   input = ArrowFile_d74,
#   k = 10, 
#   knnMethod = "UMAP", 
#   LSIMethod = 1
# )
# 
# proj_d74 <- ArchRProject(
#   ArrowFiles = ArrowFile_d74, 
#   outputDirectory = "d74",
#   copyArrows = TRUE
# )
# 
# proj_d74 <- filterDoublets(ArchRProj = proj_d74)
# 
# proj_d74 <- saveArchRProject(ArchRProj = proj_d74)

##############################################################################

# inputFile_d78 <- c("d78"="GSE183684_RAW/GSM5567520_d78_fragments.tsv.gz")
# 
# addArchRThreads(1)
# 
# ArrowFile_d78 <- createArrowFiles(
#   inputFiles = inputFile_d78,
#   sampleNames = names(inputFile_d78),
#   minTSS = 0,
#   minFrags = 1000, 
#   addTileMat = TRUE,
#   addGeneScoreMat = TRUE
# )
# 
# doubScores <- addDoubletScores(
#   input = ArrowFile_d78,
#   k = 10, 
#   knnMethod = "UMAP", 
#   LSIMethod = 1
# )
# 
# proj_d78 <- ArchRProject(
#   ArrowFiles = ArrowFile_d78, 
#   outputDirectory = "d78",
#   copyArrows = TRUE
# )
# 
# proj_d78 <- filterDoublets(ArchRProj = proj_d78)
# 
# proj_d78 <- saveArchRProject(ArchRProj = proj_d78)

##############################################################################

# inputFile_d113 <- c("d113"="GSE183684_RAW/GSM5567521_d113_fragments.tsv.gz")
# 
# addArchRThreads(1)
# 
# ArrowFile_d113 <- createArrowFiles(
#   inputFiles = inputFile_d113,
#   sampleNames = names(inputFile_d113),
#   minTSS = 0,
#   minFrags = 1000, 
#   addTileMat = TRUE,
#   addGeneScoreMat = TRUE
# )
# 
# doubScores <- addDoubletScores(
#   input = ArrowFile_d113,
#   k = 10, 
#   knnMethod = "UMAP", 
#   LSIMethod = 1
# )
# 
# proj_d113 <- ArchRProject(
#   ArrowFiles = ArrowFile_d113, 
#   outputDirectory = "d113",
#   copyArrows = TRUE
# )
# 
# proj_d113 <- filterDoublets(ArchRProj = proj_d113)
# 
# proj_d113 <- saveArchRProject(ArchRProj = proj_d113)

##############################################################################

# inputFile_d132 <- c("d132"="GSE183684_RAW/GSM5567522_d132_fragments.tsv.gz")
# 
# addArchRThreads(1)
# 
# ArrowFile_d132 <- createArrowFiles(
#   inputFiles = inputFile_d132,
#   sampleNames = names(inputFile_d132),
#   minTSS = 0,
#   minFrags = 1000, 
#   addTileMat = TRUE,
#   addGeneScoreMat = TRUE
# )
# 
# doubScores <- addDoubletScores(
#   input = ArrowFile_d132,
#   k = 10, 
#   knnMethod = "UMAP", 
#   LSIMethod = 1
# )
# 
# proj_d132 <- ArchRProject(
#   ArrowFiles = ArrowFile_d132, 
#   outputDirectory = "d132",
#   copyArrows = TRUE
# )
# 
# proj_d132 <- filterDoublets(ArchRProj = proj_d132)
# 
# proj_d132 <- saveArchRProject(ArchRProj = proj_d132)

##############################################################################



##############################################################################

# sessionInfo()
# 
# R version 4.3.1 (2023-06-16)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Sonoma 14.1.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
# 
# Random number generation:
#   RNG:     L'Ecuyer-CMRG 
#  Normal:  Inversion 
#  Sample:  Rejection 
#  
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: Asia/Tokyo
# tzcode source: internal
# 
# attached base packages:
# [1] stats4    grid      stats     graphics  grDevices utils     datasets  methods  
# [9] base     
# 
# other attached packages:
#  [1] nabor_0.5.0                       Rsamtools_2.16.0                 
#  [3] BSgenome.Hsapiens.UCSC.hg38_1.4.5 BSgenome_1.68.0                  
#  [5] rtracklayer_1.60.1                Biostrings_2.68.1                
#  [7] XVector_0.40.0                    rhdf5_2.44.0                     
#  [9] SummarizedExperiment_1.30.2       Biobase_2.60.0                   
# [11] MatrixGenerics_1.12.3             Rcpp_1.0.11                      
# [13] Matrix_1.6-1.1                    GenomicRanges_1.52.1             
# [15] GenomeInfoDb_1.36.4               IRanges_2.34.1                   
# [17] S4Vectors_0.38.2                  BiocGenerics_0.46.0              
# [19] matrixStats_1.0.0                 data.table_1.14.8                
# [21] stringr_1.5.0                     plyr_1.8.9                       
# [23] magrittr_2.0.3                    ggplot2_3.4.4                    
# [25] gtable_0.3.4                      gtools_3.9.4                     
# [27] gridExtra_2.3                     ArchR_1.0.2                      
# [29] Signac_1.11.0                     Seurat_4.9.9.9060                
# [31] SeuratObject_4.9.9.9091           sp_2.1-1                         
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.21         splines_4.3.1            later_1.3.1             
#   [4] BiocIO_1.10.0            bitops_1.0-7             tibble_3.2.1            
#   [7] polyclip_1.10-6          XML_3.99-0.14            fastDummies_1.7.3       
#  [10] lifecycle_1.0.4          globals_0.16.2           lattice_0.21-9          
#  [13] MASS_7.3-60              plotly_4.10.3            yaml_2.3.7              
#  [16] httpuv_1.6.11            sctransform_0.4.1        spam_2.9-1              
#  [19] spatstat.sparse_3.0-2    reticulate_1.34.0        cowplot_1.1.1           
#  [22] pbapply_1.7-2            RColorBrewer_1.1-3       abind_1.4-5             
#  [25] zlibbioc_1.46.0          Rtsne_0.16               purrr_1.0.2             
#  [28] RCurl_1.98-1.12          GenomeInfoDbData_1.2.10  ggrepel_0.9.4           
#  [31] irlba_2.3.5.1            listenv_0.9.0            spatstat.utils_3.0-3    
#  [34] goftest_1.2-3            RSpectra_0.16-1          spatstat.random_3.1-6   
#  [37] fitdistrplus_1.1-11      parallelly_1.36.0        leiden_0.4.3            
#  [40] codetools_0.2-19         DelayedArray_0.26.7      RcppRoll_0.3.0          
#  [43] tidyselect_1.2.0         farver_2.1.1             spatstat.explore_3.2-3  
#  [46] GenomicAlignments_1.36.0 jsonlite_1.8.7           ellipsis_0.3.2          
#  [49] progressr_0.14.0         ggridges_0.5.4           survival_3.5-7          
#  [52] tools_4.3.1              ica_1.0-3                glue_1.6.2              
#  [55] dplyr_1.1.4              withr_2.5.2              fastmap_1.1.1           
#  [58] rhdf5filters_1.12.1      fansi_1.0.5              digest_0.6.33           
#  [61] R6_2.5.1                 mime_0.12                colorspace_2.1-0        
#  [64] scattermore_1.2          tensor_1.5               spatstat.data_3.0-1     
#  [67] utf8_1.2.4               tidyr_1.3.0              generics_0.1.3          
#  [70] httr_1.4.7               htmlwidgets_1.6.2        S4Arrays_1.0.6          
#  [73] uwot_0.1.16              pkgconfig_2.0.3          lmtest_0.9-40           
#  [76] htmltools_0.5.6.1        dotCall64_1.1-0          scales_1.2.1            
#  [79] png_0.1-8                rstudioapi_0.15.0        reshape2_1.4.4          
#  [82] rjson_0.2.21             nlme_3.1-163             zoo_1.8-12              
#  [85] KernSmooth_2.23-22       parallel_4.3.1           miniUI_0.1.1.1          
#  [88] restfulr_0.0.15          pillar_1.9.0             vctrs_0.6.4             
#  [91] RANN_2.6.1               promises_1.2.1           xtable_1.8-4            
#  [94] cluster_2.1.4            cli_3.6.1                compiler_4.3.1          
#  [97] rlang_1.1.2              crayon_1.5.2             future.apply_1.11.0     
# [100] labeling_0.4.3           stringi_1.7.12           viridisLite_0.4.2       
# [103] deldir_1.0-9             BiocParallel_1.34.2      munsell_0.5.0           
# [106] lazyeval_0.2.2           spatstat.geom_3.2-7      RcppHNSW_0.5.0          
# [109] patchwork_1.1.3          future_1.33.0            Rhdf5lib_1.22.1         
# [112] shiny_1.7.5.1            ROCR_1.0-11              igraph_1.5.1            
# [115] fastmatch_1.1-4   
