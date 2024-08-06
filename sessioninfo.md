R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=English_United Kingdom.utf8  LC_CTYPE=English_United Kingdom.utf8    LC_MONETARY=English_United Kingdom.utf8 LC_NUMERIC=C                           
[5] LC_TIME=English_United Kingdom.utf8    

time zone: Europe/London
tzcode source: internal

attached base packages:
[1] stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] slingshot_2.12.0            TrajectoryUtils_1.12.0      SingleCellExperiment_1.26.0 SummarizedExperiment_1.34.0 Biobase_2.64.0              GenomicRanges_1.56.1       
 [7] GenomeInfoDb_1.40.1         IRanges_2.38.0              S4Vectors_0.42.0            BiocGenerics_0.50.0         MatrixGenerics_1.16.0       matrixStats_1.3.0          
[13] princurve_2.1.6             ggridges_0.5.6              colorRamp2_0.1.0            ComplexHeatmap_2.20.0       readxl_1.4.3                EnhancedVolcano_1.22.0     
[19] ggrepel_0.9.5               reshape2_1.4.4              ggpubr_0.6.0                sqldf_0.4-11                RSQLite_2.3.7               gsubfn_0.7                 
[25] proto_1.0.0                 stringr_1.5.1               ggplot2_3.5.1               dplyr_1.1.4                 cowplot_1.1.3               Seurat_5.1.0               
[31] SeuratObject_5.0.2          sp_2.1-4                   

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22          splines_4.4.1             later_1.3.2               R.oo_1.26.0               tibble_3.2.1              cellranger_1.1.0         
  [7] polyclip_1.10-6           fastDummies_1.7.3         lifecycle_1.0.4           tcltk_4.4.1               rstatix_0.7.2             doParallel_1.0.17        
 [13] globals_0.16.3            lattice_0.22-6            MASS_7.3-61               backports_1.5.0           magrittr_2.0.3            plotly_4.10.4            
 [19] rmarkdown_2.27            yaml_2.3.9                httpuv_1.6.15             glmGamPoi_1.16.0          sctransform_0.4.1         spam_2.10-0              
 [25] spatstat.sparse_3.1-0     reticulate_1.38.0         pbapply_1.7-2             DBI_1.2.3                 RColorBrewer_1.1-3        zlibbioc_1.50.0          
 [31] abind_1.4-5               Rtsne_0.17                R.utils_2.12.3            purrr_1.0.2               GenomeInfoDbData_1.2.12   circlize_0.4.16          
 [37] irlba_2.3.5.1             listenv_0.9.1             spatstat.utils_3.0-5      goftest_1.2-3             RSpectra_0.16-1           spatstat.random_3.2-3    
 [43] fitdistrplus_1.2-1        parallelly_1.37.1         DelayedMatrixStats_1.26.0 DelayedArray_0.30.1       leiden_0.4.3.1            codetools_0.2-20         
 [49] tidyselect_1.2.1          shape_1.4.6.1             farver_2.1.2              UCSC.utils_1.0.0          spatstat.explore_3.2-7    jsonlite_1.8.8           
 [55] GetoptLong_1.0.5          progressr_0.14.0          survival_3.7-0            iterators_1.0.14          foreach_1.5.2             tools_4.4.1              
 [61] chron_2.3-61              ica_1.0-3                 Rcpp_1.0.12               glue_1.7.0                SparseArray_1.4.8         gridExtra_2.3            
 [67] xfun_0.45                 withr_3.0.0               fastmap_1.2.0             fansi_1.0.6               digest_0.6.36             R6_2.5.1                 
 [73] mime_0.12                 colorspace_2.1-0          scattermore_1.2           tensor_1.5                spatstat.data_3.1-2       R.methodsS3_1.8.2        
 [79] utf8_1.2.4                tidyr_1.3.1               generics_0.1.3            data.table_1.15.4         S4Arrays_1.4.1            httr_1.4.7               
 [85] htmlwidgets_1.6.4         uwot_0.2.2                pkgconfig_2.0.3           gtable_0.3.5              blob_1.2.4                lmtest_0.9-40            
 [91] XVector_0.44.0            htmltools_0.5.8.1         carData_3.0-5             dotCall64_1.1-1           clue_0.3-65               scales_1.3.0             
 [97] png_0.1-8                 knitr_1.48                rstudioapi_0.16.0         rjson_0.2.21              nlme_3.1-165              cachem_1.1.0             
[103] zoo_1.8-12                GlobalOptions_0.1.2       KernSmooth_2.23-24        vipor_0.4.7               parallel_4.4.1            miniUI_0.1.1.1           
[109] ggrastr_1.0.2             pillar_1.9.0              vctrs_0.6.5               RANN_2.6.1                promises_1.3.0            car_3.1-2                
[115] xtable_1.8-4              cluster_2.1.6             beeswarm_0.4.0            evaluate_0.24.0           cli_3.6.3                 compiler_4.4.1           
[121] rlang_1.1.4               crayon_1.5.3              future.apply_1.11.2       ggsignif_0.6.4            labeling_0.4.3            ggbeeswarm_0.7.2         
[127] plyr_1.8.9                stringi_1.8.4             viridisLite_0.4.2         deldir_2.0-4              munsell_0.5.1             lazyeval_0.2.2           
[133] spatstat.geom_3.2-9       Matrix_1.7-0              RcppHNSW_0.6.0            patchwork_1.2.0           sparseMatrixStats_1.16.0  bit64_4.0.5              
[139] future_1.33.2             shiny_1.8.1.1             ROCR_1.0-11               igraph_2.0.3              broom_1.0.6               memoise_2.0.1            
[145] bit_4.0.5   