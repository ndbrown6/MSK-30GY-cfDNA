# Intra-treatment ctDNA and imaging dynamics identify candidates for adaptive therapy in HPV-related oropharyngeal cancer   

## Diplas *et al.*
### https://someurl
&nbsp;
&nbsp;
&nbsp;

![Front page](https://github.com/ndbrown6/MSK-30GY-cfDNA/blob/main/etc/splash.png)

#### Clone repository
```
git clone https://github.com/ndbrown6/MSK-30GY-cfDNA.git
cd MSK-30GY-cfDNA/R
```

#### Render markdowns in R
```
library("rmarkdown")
rmarkdown::render(input = "Figure_1__oncoplot_tcga_hnsc.Rmd", output_dir = "../res/")
```

#### R session info
```
R version 4.5.2 (2025-10-31)
Platform: aarch64-apple-darwin20.0.0
Running under: macOS Tahoe 26.5.1

Matrix products: default
BLAS:   libblas.3.9.0.dylib 
LAPACK: liblapack.3.9.0.dylib  LAPACK version 3.9.0

locale:
[1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: America/New_York
tzcode source: system (macOS)

attached base packages:
 [1] stats4    parallel  grid      stats     graphics  grDevices utils    
 [8] datasets  methods   base     

other attached packages:
 [1] survminer_0.5.2             survival_3.8-6             
 [3] randomForest_4.7-1.2        ggparty_1.0.0.1            
 [5] partykit_1.2-27             libcoin_1.0-11             
 [7] party_1.3-19                strucchange_1.5-4          
 [9] sandwich_3.1-1              zoo_1.8-15                 
[11] modeltools_0.2-24           mvtnorm_1.3-3              
[13] TCGAmutations_0.1.0         e1071_1.7-17               
[15] vip_0.4.5                   ipred_0.9-15               
[17] plyr_1.8.9                  rpart.plot_3.1.4           
[19] ggdendro_0.2.0              Rtsne_0.17                 
[21] sensemakr_0.1.6             limma_3.66.0               
[23] genefilter_1.92.0           fgsea_1.36.2               
[25] GSVAdata_1.46.0             SummarizedExperiment_1.40.0
[27] GenomicRanges_1.62.1        Seqinfo_1.0.0              
[29] MatrixGenerics_1.22.0       matrixStats_1.5.0          
[31] hgu95a.db_3.13.0            GSEABase_1.72.0            
[33] graph_1.88.1                annotate_1.88.0            
[35] XML_3.99-0.22               GSVA_2.4.7                 
[37] org.Hs.eg.db_3.22.0         AnnotationDbi_1.72.0       
[39] IRanges_2.44.0              S4Vectors_0.48.0           
[41] Biobase_2.70.0              BiocGenerics_0.56.0        
[43] generics_0.1.4              corrplot_0.95              
[45] ggbeeswarm_0.7.3            ellipse_0.5.0              
[47] PMCMRplus_1.9.12            klaR_1.7-4                 
[49] rpart_4.1.24                pROC_1.19.0.1              
[51] RColorBrewer_1.1-3          maftools_2.12.05           
[53] doMC_1.3.8                  iterators_1.0.14           
[55] foreach_1.5.2               preseqR_4.0.0              
[57] drc_3.0-1                   MASS_7.3-65                
[59] gridExtra_2.3               ComplexHeatmap_2.26.1      
[61] superheat_0.1.0             viridis_0.6.5              
[63] viridisLite_0.4.3           cowplot_1.2.0              
[65] gghalves_0.1.4              ggridges_0.5.7             
[67] ggpubr_0.6.3                ggord_1.1.8                
[69] ggforce_0.5.0               ggrepel_0.9.6              
[71] ggsignif_0.6.4              pander_0.6.6               
[73] fuzzyjoin_0.1.8             tidyr_1.3.2                
[75] magrittr_2.0.4              reshape2_1.4.5             
[77] dplyr_1.2.0                 data.table_1.17.8          
[79] readr_2.2.0                 ggplot2_4.0.2              
[81] gdata_3.0.1                 rmarkdown_2.30             

loaded via a namespace (and not attached):
  [1] httr_1.4.8                  doParallel_1.0.17          
  [3] tools_4.5.2                 backports_1.5.0            
  [5] R6_2.6.1                    HDF5Array_1.38.0           
  [7] questionr_0.8.2             rhdf5filters_1.22.0        
  [9] GetoptLong_1.1.0            withr_3.0.2                
 [11] cli_3.6.5                   labeling_0.4.3             
 [13] sass_0.4.10                 BWStest_0.2.3              
 [15] S7_0.2.1                    proxy_0.4-29               
 [17] dichromat_2.0-0.1           parallelly_1.46.1          
 [19] labelled_2.16.0             plotrix_3.8-14             
 [21] rstudioapi_0.18.0           RSQLite_2.4.6              
 [23] shape_1.4.6.1               combinat_0.0-8             
 [25] vroom_1.7.0                 gtools_3.9.5               
 [27] car_3.1-5                   Matrix_1.7-4               
 [29] abind_1.4-8                 lifecycle_1.0.5            
 [31] multcomp_1.4-30             yaml_2.3.12                
 [33] inum_1.0-5                  carData_3.0-6              
 [35] rhdf5_2.54.1                SparseArray_1.10.8         
 [37] blob_1.3.0                  promises_1.5.0             
 [39] crayon_1.5.3                miniUI_0.1.2               
 [41] lattice_0.22-9              haven_2.5.5                
 [43] beachmat_2.26.0             KEGGREST_1.50.0            
 [45] magick_2.9.1                pillar_1.11.1              
 [47] knitr_1.51                  rjson_0.2.23               
 [49] future.apply_1.20.2         kSamples_1.2-12            
 [51] codetools_0.2-20            fastmatch_1.1-8            
 [53] glue_1.8.0                  memuse_4.2-3               
 [55] vctrs_0.7.1                 png_0.1-8                  
 [57] gtable_0.3.6                cachem_1.1.0               
 [59] xfun_0.56                   prodlim_2025.04.28         
 [61] S4Arrays_1.10.1             mime_0.13                  
 [63] SingleCellExperiment_1.32.0 lava_1.8.2                 
 [65] statmod_1.5.1               gmp_0.7-5.1                
 [67] TH.data_1.1-5               bit64_4.6.0-1              
 [69] bslib_0.10.0                irlba_2.3.7                
 [71] vipor_0.4.7                 otel_0.2.0                 
 [73] colorspace_2.1-2            DBI_1.3.0                  
 [75] nnet_7.3-20                 DNAcopy_1.84.0             
 [77] tidyselect_1.2.1            bit_4.6.0                  
 [79] compiler_4.5.2              h5mread_1.2.1              
 [81] DelayedArray_0.36.0         checkmate_2.3.3            
 [83] scales_1.4.0                multcompView_0.1-11        
 [85] stringr_1.6.0               SpatialExperiment_1.20.0   
 [87] digest_0.6.39               XVector_0.50.0             
 [89] htmltools_0.5.9             pkgconfig_2.0.3            
 [91] sparseMatrixStats_1.22.0    highr_0.12                 
 [93] fastmap_1.2.0               rlang_1.1.7                
 [95] GlobalOptions_0.1.3         shiny_1.13.0               
 [97] SuppDists_1.1-9.9           farver_2.1.2               
 [99] jquerylib_0.1.4             jsonlite_2.0.0             
[101] BiocParallel_1.44.0         BiocSingular_1.26.1        
[103] polynom_1.4-1               Formula_1.2-5              
[105] Rhdf5lib_1.32.0             Rcpp_1.1.1                 
[107] stringi_1.8.7               listenv_0.10.0             
[109] forcats_1.0.1               Biostrings_2.78.0          
[111] splines_4.5.2               hms_1.1.4                  
[113] circlize_0.4.17             ScaledMatrix_1.18.0        
[115] evaluate_1.0.5              tzdb_0.5.0                 
[117] tweenr_2.0.3                httpuv_1.6.16              
[119] purrr_1.2.1                 polyclip_1.10-7            
[121] future_1.69.0               clue_0.3-67                
[123] coin_1.4-3                  rsvd_1.0.5                 
[125] broom_1.0.12                xtable_1.8-8               
[127] Rmpfr_1.1-1                 rstatix_0.7.3              
[129] later_1.4.8                 class_7.3-23               
[131] tibble_3.3.1                memoise_2.0.1              
[133] beeswarm_0.4.0              cluster_2.1.8.2            
[135] globals_0.19.0  
```