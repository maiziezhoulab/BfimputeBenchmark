# BfimputePaper
This repository contains the functions that are used to generate the results of
paper Bfimpute。

The Bfimpute source code can be found [here](https://github.com/maiziezhoulab/BfimputePaper).

Questions about the paper or scripts can be sent to wenzihang0506@gmail.com.


## Session Information
<details>
<summary><code>utils::sessionInfo()</code></summary>
<pre>
R version 3.6.1 (2019-07-05)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /share/software/user/open/openblas/0.2.19/lib/libopenblasp-r0.2.19.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] splines   parallel  stats4    stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] org.Hs.eg.db_3.10.0         AnnotationDbi_1.48.0       
 [3] clusterProfiler_3.14.3      WGCNA_1.70-3               
 [5] fastcluster_1.2.3           dynamicTreeCut_1.63-1      
 [7] monocle_2.14.0              DDRTree_0.1.5              
 [9] irlba_2.3.3                 VGAM_1.1-5                 
[11] Matrix_1.2-17               TSCAN_1.7.0                
[13] ggbeeswarm_0.6.0            ggpubr_0.4.0               
[15] pheatmap_1.0.12             DESeq2_1.26.0              
[17] cluster_2.1.0               Spectrum_1.1               
[19] ggthemes_4.2.4              gridExtra_2.3              
[21] cowplot_1.1.1               scater_1.14.6              
[23] ggplot2_3.3.5               SingleCellExperiment_1.8.0 
[25] SummarizedExperiment_1.16.1 DelayedArray_0.12.3        
[27] BiocParallel_1.20.1         matrixStats_0.60.0         
[29] Biobase_2.46.0              GenomicRanges_1.38.0       
[31] GenomeInfoDb_1.22.1         IRanges_2.20.2             
[33] S4Vectors_0.24.4            BiocGenerics_0.32.0        

loaded via a namespace (and not attached):
  [1] ClusterR_1.2.5           tidyr_1.1.3              bit64_4.0.5             
  [4] knitr_1.33               data.table_1.14.0        rpart_4.1-15            
  [7] RCurl_1.98-1.3           doParallel_1.0.16        generics_0.1.0          
 [10] preprocessCore_1.48.0    RSQLite_2.2.7            RANN_2.6.1              
 [13] europepmc_0.4            combinat_0.0-8           bit_4.0.4               
 [16] enrichplot_1.6.1         xml2_1.3.2.9001          httpuv_1.6.1            
 [19] assertthat_0.2.1         viridis_0.6.1            xfun_0.25               
 [22] hms_1.1.0                promises_1.2.0.1         fansi_0.5.0             
 [25] progress_1.2.2           caTools_1.18.2           readxl_1.3.1            
 [28] igraph_1.2.6             DBI_1.1.1                geneplotter_1.64.0      
 [31] htmlwidgets_1.5.3        sparsesvd_0.2            purrr_0.3.4             
 [34] ellipsis_0.3.2           dplyr_1.0.7              backports_1.2.1         
 [37] annotate_1.64.0          vctrs_0.3.8              abind_1.4-5             
 [40] cachem_1.0.5             withr_2.4.2              ggforce_0.3.3           
 [43] triebeard_0.3.0          checkmate_2.0.0          prettyunits_1.1.1       
 [46] mclust_5.4.7             DOSE_3.12.0              crayon_1.4.1            
 [49] genefilter_1.68.0        pkgconfig_2.0.3          slam_0.1-48             
 [52] tweenr_1.0.2             nlme_3.1-140             vipor_0.4.5             
 [55] nnet_7.3-12              rlang_0.4.11             diptest_0.76-0          
 [58] lifecycle_1.0.0          rsvd_1.0.3               cellranger_1.1.0        
 [61] polyclip_1.10-0          urltools_1.7.3           carData_3.0-4           
 [64] base64enc_0.1-3          beeswarm_0.4.0           ggridges_0.5.3          
 [67] png_0.1-7                viridisLite_0.4.0        bitops_1.0-7            
 [70] KernSmooth_2.23-15       blob_1.2.2               DelayedMatrixStats_1.8.0
 [73] stringr_1.4.0            qvalue_2.18.0            jpeg_0.1-9              
 [76] rstatix_0.7.0            gridGraphics_0.5-1       ggsignif_0.6.2          
 [79] scales_1.1.1             memoise_2.0.0            magrittr_2.0.1          
 [82] plyr_1.8.6               gplots_3.1.1             zlibbioc_1.32.0         
 [85] compiler_3.6.1           HSMMSingleCell_1.6.0     RColorBrewer_1.1-2      
 [88] XVector_0.26.0           htmlTable_2.2.1          Formula_1.2-4           
 [91] MASS_7.3-51.4            mgcv_1.8-28              tidyselect_1.1.1        
 [94] stringi_1.7.3            forcats_0.5.1            densityClust_0.3        
 [97] GOSemSim_2.12.1          BiocSingular_1.2.2       locfit_1.5-9.4          
[100] latticeExtra_0.6-29      ggrepel_0.9.1            grid_3.6.1              
[103] fastmatch_1.1-3          tools_3.6.1              rio_0.5.27              
[106] rstudioapi_0.13          foreach_1.5.1            foreign_0.8-71          
[109] farver_2.1.0             Rtsne_0.15               ggraph_2.0.5            
[112] RcppZiggurat_0.1.6       digest_0.6.27            rvcheck_0.1.8           
[115] BiocManager_1.30.16      FNN_1.1.3                shiny_1.6.0             
[118] qlcMatrix_0.9.7          Rcpp_1.0.7               car_3.0-11              
[121] broom_0.7.9              later_1.2.0              httr_1.4.2              
[124] colorspace_2.0-2         XML_3.99-0.3             graphlayouts_0.7.1      
[127] ggplotify_0.0.8          xtable_1.8-4             gmp_0.6-1               
[130] jsonlite_1.7.2           tidygraph_1.2.0          Rfast_2.0.3             
[133] R6_2.5.0                 Hmisc_4.5-0              pillar_1.6.2            
[136] htmltools_0.5.1.1        mime_0.11                glue_1.4.2              
[139] fastmap_1.1.0            BiocNeighbors_1.4.2      codetools_0.2-16        
[142] fgsea_1.12.0             utf8_1.2.2               lattice_0.20-38         
[145] tibble_3.1.3             curl_4.3.2               gtools_3.9.2            
[148] zip_2.2.0                GO.db_3.10.0             openxlsx_4.2.4          
[151] survival_3.2-12          limma_3.42.2             docopt_0.7.1            
[154] fastICA_1.2-2            munsell_0.5.0            DO.db_2.9               
[157] GenomeInfoDbData_1.2.2   iterators_1.0.13         impute_1.60.0           
[160] haven_2.4.3              reshape2_1.4.4           gtable_0.3.0            
</pre>
</details>
