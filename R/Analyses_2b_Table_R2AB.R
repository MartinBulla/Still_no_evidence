# ===========================================================================
# Supporting information for "Bulla, Valcu & Kempenaers. (2021) "Still no 
# evidence for disruption of global patterns of nest predation in shorebirds" 
# Contributor: Martin Bulla
# 📍 script runs relative to the project's root directory, uses Kubelka et 
# al's data to generates Table R1AB within part
# (2) Testing the key interaction effect: does the temporal change in 
# predation vary across the globe? 
# ===========================================================================

# load tools and data
  source('R/Constants_Functions.R')

  d <- read.csv( "Data/DATApopulations.csv", h = T, sep=";")   
  d$site = paste(d$Latitude,d$Longitude)
  d$lat_abs = abs(d$Latitude) # abs latitude
  d$ln_N_nests = log(d$N_nests)
  d$hemisphere =as.factor(ifelse(d$Latitude > 0, "Northern", "Southern"))

# DPR - 2 way as in Kubelka et al 2018a
  md2= lmer(log(DPR) ~ log(N_nests) + hemisphere+scale(mean_year)+scale(abs(Latitude)) + hemisphere:scale(mean_year)+hemisphere:scale(abs(Latitude)) +scale(mean_year):scale(abs(Latitude))+(1|site)+(1|species),  data = d) 
  summary(glht(md2))

# TPR - 2 way as in Kubelka et al 2018a
  mt2= lmer(TPR ~ log(N_nests) + hemisphere+scale(mean_year)+scale(abs(Latitude)) + hemisphere:scale(mean_year)+hemisphere:scale(abs(Latitude)) +scale(mean_year):scale(abs(Latitude))+(1|site)+(1|species),  data = d) 
  summary(glht(mt2))

# DPR - 3 way as in Kubelka et al 2018b script
  md3= lmer(log(DPR) ~ log( N_nests) + hemisphere*scale(mean_year)*scale(abs(Latitude))  +(1|site)+(1|species),  data = d) 
  summary(glht(md3))

# TPR - 3 way as in Kubelka et al 2018b script
  mt3= lmer(TPR ~ log( N_nests) + hemisphere*scale(mean_year)*scale(abs(Latitude))  +(1|site)+(1|species),  data = d) 
  summary(glht(mt3))

# export results
  o1 = m_out(name = "dpr 2-way", model = md2, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
  o2 = m_out(name = "tpr 2-way", model = mt2, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
  o3 = m_out(name = "dpr 3-way", model = md3, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
  o4 = m_out(name = "tpr 3-way", model = mt3, round_ = 3, nsim = 5000, aic = FALSE)[1:6]

  o1$p_value = NA
  o1$p_value[1:8] = round(summary(glht(md2))$test$pvalues,2)
  o2$p_value = NA
  o2$p_value[1:8] = round(summary(glht(mt2))$test$pvalues,2)
  o3$p_value = NA
  o3$p_value[1:9] = round(summary(glht(md3))$test$pvalues,2)
  o4$p_value = NA
  o4$p_value[1:9] = round(summary(glht(mt3))$test$pvalues,2)

  
  write_xlsx(rbind(o1,o2,o3,o4), 'Outputs/Table_R2AB.xlsx')

# sessionInfo()
  # R version 4.0.2 (2020-06-22)
  # Platform: x86_64-apple-darwin17.0 (64-bit)
  # Running under: macOS Mojave 10.14.6
  # 
  # Matrix products: default
  # BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
  # LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
  # 
  # locale:
  # [1] C/UTF-8/C/C/C/C
  # 
  # attached base packages:
  # [1] stats     graphics  grDevices utils     datasets  methods   base     
  # 
  # other attached packages:
  #  [1] writexl_1.3.1      performance_0.4.8  RColorBrewer_1.1-2 multcomp_1.4-13    TH.data_1.0-10    
  #  [6] survival_3.1-12    mvtnorm_1.1-1      arm_1.11-2         lme4_1.1-23        Matrix_1.2-18     
  # [11] MASS_7.3-51.6     
  # 
  # loaded via a namespace (and not attached):
  #  [1] Rcpp_1.0.5          lattice_0.20-41     png_0.1-7           zoo_1.8-8          
  #  [5] digest_0.6.25       R6_2.4.1            backports_1.1.8     acepack_1.4.1      
  #  [9] coda_0.19-3         ggplot2_3.3.2       pillar_1.4.6        rlang_0.4.7        
  # [13] rstudioapi_0.11     minqa_1.2.4         data.table_1.13.0   nloptr_1.2.2.2     
  # [17] rpart_4.1-15        checkmate_2.0.0     splines_4.0.2       statmod_1.4.34     
  # [21] stringr_1.4.0       foreign_0.8-80      htmlwidgets_1.5.1   munsell_0.5.0      
  # [25] compiler_4.0.2      xfun_0.16           pkgconfig_2.0.3     base64enc_0.1-3    
  # [29] htmltools_0.5.0     nnet_7.3-14         insight_0.9.0       tidyselect_1.1.0   
  # [33] tibble_3.0.3        gridExtra_2.3       htmlTable_2.0.1     Hmisc_4.4-0        
  # [37] codetools_0.2-16    crayon_1.3.4        dplyr_1.0.1         grid_4.0.2         
  # [41] nlme_3.1-148        gtable_0.3.0        lifecycle_0.2.0     magrittr_1.5       
  # [45] bayestestR_0.7.2    scales_1.1.1        stringi_1.5.3       latticeExtra_0.6-29
  # [49] ellipsis_0.3.1      generics_0.0.2      vctrs_0.3.2         boot_1.3-25        
  # [53] sandwich_2.5-1      Formula_1.2-3       tools_4.0.2         glue_1.4.2         
  # [57] purrr_0.3.4         jpeg_0.1-8.1        abind_1.4-5         colorspace_1.4-1   
  # [61] cluster_2.1.0       knitr_1.29  