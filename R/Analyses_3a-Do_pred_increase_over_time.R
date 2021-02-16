# ===========================================================================
# Supporting information for "Bulla, Valcu & Kempenaers. (2021) "Still no 
# evidence for disruption of global patterns of nest predation in shorebirds" 
# Contributor: Martin Bulla
# 📍 script runs relative to the project's root directory, uses Kubelka et 
# al's data to generates outputs for part
# (3) Do predation rates increase over time? 
# ===========================================================================

# load packages, constants and data
    rm( list = ls() )
    source('R/Constants_Functions.R')
    d = read.csv("Data/DATApopulations.csv", h = T, sep=";",stringsAsFactors = FALSE)
    d$site = paste(d$Latitude,d$Longitude) # define site
    dn = d[d$DPR_tran == 'NO',] # only high quality data
  
# Kubelka et al.'s (2018) Table S3a model	
     mk = lm(log(DPR) ~ log( N_nests) + mean_year,  data = d)
     summary(mk)
# Kubelka et al.'s (2018) Table S3a model CONTROLLED FOR STUDY SITE - effect reduced by 30%
     mkc = lmer(log(DPR) ~ log( N_nests) + mean_year+(1|site),  data = d) 
     summary(mkc)

# Kubelka et al.'s (2018) Table S3d model on non-transformed data 
     mkn = lm(log(DPR) ~ log( N_nests) +mean_year,  data = dn)  
     summary(mkn)
# Kubelka et al.'s (2018) Table S3d model on non-transformed data CONTROLLED FOR STUDY SITE - effect reduced by 30%
     mknc = lmer(log(DPR) ~ log( N_nests) +mean_year   +(1|site),  data = dn)  
     summary(mknc)
     summary(glht(mknc))
     m_out(name = "kmnc", model = mknc, round_ = 3, nsim = 5000, aic = FALSE) 

# sessionInfo()
  #R version 4.0.2 (2020-06-22)
  #Platform: x86_64-apple-darwin17.0 (64-bit)
  #Running under: macOS Mojave 10.14.6
  #
  #Matrix products: default
  #BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
  #LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
  #
  #locale:
  #[1] C/UTF-8/C/C/C/C
  #
  #attached base packages:
  #[1] stats     graphics  grDevices utils     datasets  methods   base     
  #
  #other attached packages:
  # [1] ggplot2_3.3.2      writexl_1.3.1      RColorBrewer_1.1-2 performance_0.4.8  multcomp_1.4-13   
  # [6] TH.data_1.0-10     survival_3.1-12    mvtnorm_1.1-1      data.table_1.13.0  arm_1.11-2        
  #[11] lme4_1.1-23        Matrix_1.2-18      MASS_7.3-51.6     
  #
  #loaded via a namespace (and not attached):
  # [1] Rcpp_1.0.5          lattice_0.20-41     png_0.1-7           zoo_1.8-8           digest_0.6.25      
  # [6] R6_2.4.1            backports_1.1.8     acepack_1.4.1       coda_0.19-3         pillar_1.4.6       
  #[11] rlang_0.4.7         rstudioapi_0.11     minqa_1.2.4         nloptr_1.2.2.2      rpart_4.1-15       
  #[16] checkmate_2.0.0     splines_4.0.2       statmod_1.4.34      stringr_1.4.0       foreign_0.8-80     
  #[21] htmlwidgets_1.5.1   munsell_0.5.0       compiler_4.0.2      xfun_0.16           pkgconfig_2.0.3    
  #[26] base64enc_0.1-3     htmltools_0.5.0     nnet_7.3-14         insight_0.9.0       tidyselect_1.1.0   
  #[31] tibble_3.0.3        gridExtra_2.3       htmlTable_2.0.1     Hmisc_4.4-0         codetools_0.2-16   
  #[36] withr_2.2.0         crayon_1.3.4        dplyr_1.0.1         grid_4.0.2          nlme_3.1-148       
  #[41] gtable_0.3.0        lifecycle_0.2.0     magrittr_1.5        bayestestR_0.7.2    scales_1.1.1       
  #[46] stringi_1.5.3       latticeExtra_0.6-29 ellipsis_0.3.1      generics_0.0.2      vctrs_0.3.2        
  #[51] boot_1.3-25         sandwich_2.5-1      Formula_1.2-3       tools_4.0.2         glue_1.4.2         
  #[56] purrr_0.3.4         jpeg_0.1-8.1        abind_1.4-5         colorspace_1.4-1    cluster_2.1.0      
  #[61] knitr_1.29  