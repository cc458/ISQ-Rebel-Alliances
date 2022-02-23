# Replication Data for: "Do Birds of a Feather Flock Together? Rebel Constituencies and Civil War Alliances""


## Please use the following citation: 

Laia Balcells, Chong Chen, and Costantino Pischedda , “Do Birds of a Feather Flock Together? Rebel Constituencies and Civil War Alliances.” *International Studies Quarterly*, Volume 66, Issue 1, March 2022, sqab095, https://doi.org/10.1093/isq/sqab095

 * For questions about the data, please contact the authors: Laia Balcells (laia.balcells@georgetown.edu), Chong Chen (chongchen@tsinghua.edu.cn) and Costantino Pischedda(cpischedda@miami.edu).
 * Please note that there are large files that are uploaded with `git lfs` in this repository. To make sure you get the right size, we suggest that you use the `git lfs pull` or `git lfs clone` command to download the entire repository.



# Replication Instruction: 

- Please check the `ReadMe` file first. There are four folders within the repository.


- **1. Data**:  The `data` folder includes data files for all models
   + `glm_df4.RData`: Main data used for all of our models except for the AME Models
   + `ame`: A folder for AME models and the model outputs (Figure 5 and Figure A11)
		+ ``AllyList.RData``: Network data  for rebel alliances 
		+ ``Formal_AllyList.RData``: Network data  for rebel formal alliances 
		+ ``xDyadL.RData``: Dyadic-level covarites for AME
		+ ``xNodeL.RData``: Nodal-level  covarites for AME
		+ ``Fit_AME47_15.RData``: AME results for 1947-2015 (used for making Figure 5 in the main text)
		+ ``Fit_AME7909.RData``: AME results for 1947-2015 (used for making Figure A11 in the appendix)


- **2. Code**: The R code for all figures and analysis; Please run them in the following order: 

	``0_Figures 1-2.R``, ``1_Figures 3-5.R``, and ``2_Appendix_figures A1-A21.R``

	- (1) ``0_Figures 1-2.R``: R code for Figure 1-2 in the text
        - (2) ``1_Figures 3-5.R``: R code for Figure 3-5 in the text
		+ ``1_Running_AME.R``: R code for AME models (Figure 5 and Figure A11)
        - (3) “2_Appendix_figures A1-A21”: Robustness check and all appendix figures
		+ ``function_interaction_plot.R``:  Help R functions
		+ ``helpFUN.R``: help R functions; They will be loaded when running the above code
                      + ``spdur_figure_setup.R``: Help R functions

- **3. Figure**:  a folder to store all the figures in the main text

- **4. Appendix_figures**: a folder to store all the figures in the appendix



## R sessionInfo()

R version 4.0.2 (2020-06-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Mojave 10.14.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
 [1] grid      parallel  splines   stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] cowplot_1.1.1          stringi_1.5.3          ggmap_3.0.0            scales_1.1.1          
 [5] cshapes_0.6            plyr_1.8.6             maptools_1.1-1         sp_1.4-5              
 [9] directlabels_2020.6.17 ggrepel_0.9.1          ggthemes_4.2.0         ggridges_0.5.2        
[13] viridis_0.5.1          viridisLite_0.3.0      arm_1.11-2             lme4_1.1-26           
[17] Matrix_1.2-18          MASS_7.3-51.6          postregplots_0.1.0     stringr_1.4.0         
[21] gridExtra_2.3          lmtest_0.9-38          zoo_1.8-8              multiwayvcov_1.2.3    
[25] spduration_0.17.1      tidyr_1.1.2            readstata13_0.9.2      doParallel_1.0.16     
[29] iterators_1.0.12       foreach_1.5.0          reshape2_1.4.4         abind_1.4-5           
[33] texreg_1.37.5          purrr_0.3.4            VGAM_1.1-3             statnet_2019.6        
[37] tsna_0.3.1             sna_2.5                statnet.common_4.4.1   ergm.count_3.4.0      
[41] tergm_3.7.0            networkDynamic_0.10.1  ergm_3.11.0            network_1.16.0        
[45] amen_1.4.4             ggplot2_3.3.3          countrycode_1.2.0      readxl_1.3.1          
[49] readr_1.4.0            dplyr_1.0.3            haven_2.3.1           

loaded via a namespace (and not attached):
 [1] backports_1.2.1     Hmisc_4.4-1         digest_0.6.27       htmltools_0.5.0     separationplot_1.3 
 [6] fansi_0.4.2         magrittr_2.0.1      checkmate_2.0.0     rle_0.9.2           cluster_2.1.0      
[11] xts_0.12-0          sandwich_3.0-0      forecast_8.13       tseries_0.10-47     lpSolve_5.6.15     
[16] jpeg_0.1-8.1        colorspace_2.0-0    rgdal_1.5-16        xfun_0.20           crayon_1.4.1       
[21] survival_3.2-7      glue_1.4.2          gtable_0.3.0        quantmod_0.4.17     DEoptimR_1.0-8     
[26] DBI_1.1.0           Rcpp_1.0.6          xtable_1.8-4        htmlTable_2.0.1     foreign_0.8-80     
[31] Formula_1.2-3       htmlwidgets_1.5.2   httr_1.4.2          RColorBrewer_1.1-2  ellipsis_0.3.1     
[36] pkgconfig_2.0.3     farver_2.1.0        nnet_7.3-14         utf8_1.2.1          tidyselect_1.1.0   
[41] labeling_0.4.2      rlang_0.4.10        munsell_0.5.0       cellranger_1.1.0    tools_4.0.2        
[46] cli_2.3.1           generics_0.1.0      yaml_2.2.1          knitr_1.31          robustbase_0.93-6  
[51] RgoogleMaps_1.4.5.3 nlme_3.1-148        compiler_4.0.2      rstudioapi_0.13     curl_4.3           
[56] png_0.1-7           tibble_3.1.0        statmod_1.4.35      rgeos_0.5-3         forcats_0.5.1      
[61] lattice_0.20-41     nloptr_1.2.2.2      urca_1.3-0          vctrs_0.3.6         pillar_1.5.1       
[66] lifecycle_1.0.0     trust_0.1-8         corpcor_1.6.9       bitops_1.0-6        data.table_1.13.6  
[71] R6_2.5.0            latticeExtra_0.6-29 codetools_0.2-16    boot_1.3-25         assertthat_0.2.1   
[76] rjson_0.2.20        withr_2.4.1         fracdiff_1.5-1      hms_1.0.0           quadprog_1.5-8     
[81] rpart_4.1-15        timeDate_3043.102   coda_0.19-4         minqa_1.2.4         TTR_0.24.2         
[86] base64enc_0.1-3 
