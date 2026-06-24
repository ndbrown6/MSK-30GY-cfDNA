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
rmd_files = dir(path = ".", pattern = ".Rmd")
for (i in 1:length(rmd_files)) {
	rmarkdown::render(input = rmd_files[i], output_dir = "../res/")
}
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
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] fuzzyjoin_0.1.8     corrplot_0.95       klaR_1.7-4         
 [4] broom_1.0.12        TCGAmutations_0.1.0 RColorBrewer_1.1-3 
 [7] maftools_2.12.05    MASS_7.3-65         ggpubr_0.6.3       
[10] ggsignif_0.6.4      reshape2_1.4.5      magrittr_2.0.4     
[13] tidyr_1.3.2         dplyr_1.2.0         readr_2.2.0        
[16] ggplot2_4.0.2       rmarkdown_2.30     

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1    farver_2.1.2        S7_0.2.1           
 [4] fastmap_1.2.0       combinat_0.0-8      promises_1.5.0     
 [7] labelled_2.16.0     digest_0.6.39       mime_0.13          
[10] lifecycle_1.0.5     survival_3.8-6      compiler_4.5.2     
[13] rlang_1.1.7         sass_0.4.10         tools_4.5.2        
[16] yaml_2.3.12         data.table_1.17.8   knitr_1.51         
[19] labeling_0.4.3      bit_4.6.0           plyr_1.8.9         
[22] abind_1.4-8         miniUI_0.1.2        withr_3.0.2        
[25] purrr_1.2.1         BiocGenerics_0.56.0 stats4_4.5.2       
[28] grid_4.5.2          xtable_1.8-8        scales_1.4.0       
[31] dichromat_2.0-0.1   cli_3.6.5           crayon_1.5.3       
[34] generics_0.1.4      otel_0.2.0          rstudioapi_0.18.0  
[37] tzdb_0.5.0          DNAcopy_1.84.0      cachem_1.1.0       
[40] stringr_1.6.0       splines_4.5.2       parallel_4.5.2     
[43] vctrs_0.7.1         Matrix_1.7-4        jsonlite_2.0.0     
[46] carData_3.0-6       car_3.1-5           S4Vectors_0.48.0   
[49] IRanges_2.44.0      hms_1.1.4           bit64_4.6.0-1      
[52] rstatix_0.7.3       Formula_1.2-5       jquerylib_0.1.4    
[55] glue_1.8.0          cowplot_1.2.0       stringi_1.8.7      
[58] gtable_0.3.6        questionr_0.8.2     later_1.4.8        
[61] tibble_3.3.1        pillar_1.11.1       htmltools_0.5.9    
[64] R6_2.6.1            vroom_1.7.0         evaluate_1.0.5     
[67] shiny_1.13.0        lattice_0.22-9      haven_2.5.5        
[70] highr_0.12          backports_1.5.0     httpuv_1.6.16      
[73] bslib_0.10.0        Rcpp_1.1.1          gridExtra_2.3      
[76] nlme_3.1-168        mgcv_1.9-4          xfun_0.56          
[79] forcats_1.0.1       pkgconfig_2.0.3 
```