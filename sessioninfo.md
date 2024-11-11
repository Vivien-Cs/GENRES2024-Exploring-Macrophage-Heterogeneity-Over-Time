R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8   
[3] LC_MONETARY=English_United Kingdom.utf8 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.utf8    

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] clustree_0.5.1              ggraph_2.2.1                harmony_1.2.1              
 [4] Rcpp_1.0.13                 openxlsx_4.2.7.1            tidyr_1.3.1                
 [7] circlize_0.4.16             ComplexHeatmap_2.20.0       monocle3_1.3.7             
[10] viridis_0.6.5               viridisLite_0.4.2           slingshot_2.12.0           
[13] TrajectoryUtils_1.12.0      SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0
[16] Biobase_2.64.0              GenomicRanges_1.56.1        GenomeInfoDb_1.40.1        
[19] IRanges_2.38.1              S4Vectors_0.42.1            BiocGenerics_0.50.0        
[22] MatrixGenerics_1.16.0       matrixStats_1.4.1           princurve_2.1.6            
[25] readxl_1.4.3                EnhancedVolcano_1.22.0      ggrepel_0.9.6              
[28] reshape2_1.4.4              ggpubr_0.6.0                sqldf_0.4-11               
[31] RSQLite_2.3.7               gsubfn_0.7                  proto_1.0.0                
[34] stringr_1.5.1               ggplot2_3.5.1               dplyr_1.1.4                
[37] cowplot_1.1.3               Seurat_5.1.0                SeuratObject_5.0.2         
[40] sp_2.1-4                   

loaded via a namespace (and not attached):
  [1] spatstat.sparse_3.1-0     httr_1.4.7                RColorBrewer_1.1-3        doParallel_1.0.17        
  [5] tools_4.4.2               sctransform_0.4.1         backports_1.5.0           utf8_1.2.4               
  [9] R6_2.5.1                  lazyeval_0.2.2            uwot_0.2.2                GetoptLong_1.0.5         
 [13] withr_3.0.1               prettyunits_1.2.0         gridExtra_2.3             progressr_0.14.0         
 [17] cli_3.6.3                 Cairo_1.6-2               spatstat.explore_3.3-2    fastDummies_1.7.4        
 [21] labeling_0.4.3            spatstat.data_3.1-2       ggridges_0.5.6            pbapply_1.7-2            
 [25] R.utils_2.12.3            parallelly_1.38.0         limma_3.60.6              rstudioapi_0.16.0        
 [29] generics_0.1.3            shape_1.4.6.1             ica_1.0-3                 spatstat.random_3.3-2    
 [33] car_3.1-3                 zip_2.3.1                 Matrix_1.7-1              ggbeeswarm_0.7.2         
 [37] fansi_1.0.6               abind_1.4-8               R.methodsS3_1.8.2         lifecycle_1.0.4          
 [41] yaml_2.3.10               carData_3.0-5             SparseArray_1.4.8         Rtsne_0.17               
 [45] glmGamPoi_1.16.0          blob_1.2.4                promises_1.3.0            crayon_1.5.3             
 [49] miniUI_0.1.1.1            lattice_0.22-6            magick_2.8.5              pillar_1.9.0             
 [53] knitr_1.48                tcltk_4.4.2               rjson_0.2.23              boot_1.3-31              
 [57] future.apply_1.11.2       codetools_0.2-20          leiden_0.4.3.1            glue_1.8.0               
 [61] spatstat.univar_3.0-1     data.table_1.16.0         vctrs_0.6.5               png_0.1-8                
 [65] spam_2.11-0               cellranger_1.1.0          gtable_0.3.5              cachem_1.1.0             
 [69] xfun_0.48                 S4Arrays_1.4.1            mime_0.12                 tidygraph_1.3.1          
 [73] survival_3.7-0            iterators_1.0.14          statmod_1.5.0             fitdistrplus_1.2-1       
 [77] ROCR_1.0-11               nlme_3.1-166              bit64_4.5.2               progress_1.2.3           
 [81] RcppAnnoy_0.0.22          irlba_2.3.5.1             vipor_0.4.7               KernSmooth_2.23-24       
 [85] colorspace_2.1-1          DBI_1.2.3                 ggrastr_1.0.2             tidyselect_1.2.1         
 [89] bit_4.5.0                 compiler_4.4.2            chron_2.3-61              DelayedArray_0.30.1      
 [93] plotly_4.10.4             checkmate_2.3.2           scales_1.3.0              lmtest_0.9-40            
 [97] digest_0.6.36             goftest_1.2-3             presto_1.0.0              spatstat.utils_3.1-0     
[101] minqa_1.2.8               rmarkdown_2.28            XVector_0.44.0            RhpcBLASctl_0.23-42      
[105] htmltools_0.5.8.1         pkgconfig_2.0.3           lme4_1.1-35.5             sparseMatrixStats_1.16.0 
[109] fastmap_1.2.0             rlang_1.1.4               GlobalOptions_0.1.2       htmlwidgets_1.6.4        
[113] UCSC.utils_1.0.0          shiny_1.9.1               DelayedMatrixStats_1.26.0 farver_2.1.2             
[117] zoo_1.8-12                jsonlite_1.8.9            R.oo_1.26.0               magrittr_2.0.3           
[121] Formula_1.2-5             GenomeInfoDbData_1.2.12   dotCall64_1.2             patchwork_1.3.0          
[125] munsell_0.5.1             reticulate_1.39.0         stringi_1.8.4             zlibbioc_1.50.0          
[129] MASS_7.3-61               MAST_1.30.0               plyr_1.8.9                parallel_4.4.2           
[133] listenv_0.9.1             deldir_2.0-4              graphlayouts_1.2.0        splines_4.4.2            
[137] tensor_1.5                hms_1.1.3                 igraph_2.0.3              spatstat.geom_3.3-3      
[141] ggsignif_0.6.4            RcppHNSW_0.6.0            evaluate_1.0.0            tweenr_2.0.3             
[145] nloptr_2.1.1              foreach_1.5.2             httpuv_1.6.15             RANN_2.6.2               
[149] purrr_1.0.2               polyclip_1.10-7           future_1.34.0             clue_0.3-65              
[153] scattermore_1.2           ggforce_0.4.2             broom_1.0.7               xtable_1.8-4             
[157] RSpectra_0.16-2           rstatix_0.7.2             later_1.3.2               tibble_3.2.1             
[161] memoise_2.0.1             beeswarm_0.4.0            cluster_2.1.6             globals_0.16.3 