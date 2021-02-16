# ===========================================================================
# Supporting information for "Bulla, Valcu & Kempenaers. (2021) "Still no 
# evidence for disruption of global patterns of nest predation in shorebirds" 
# Contributor: Martin Bulla
# 📍 script runs relative to the project's root directory and uses Kubelka et 
# al's data and Bulla et al 2019 data preparation scripts to generate 
# Note's 6 Figure and Table
# ===========================================================================

# load packages, constants and data
    rm( list = ls() ) 
    source('R/Constants_Functions.R')
    require(ape); require(readxl)
    ax_lines = "grey60" 

    wd = here::here('Data/') # to run the Bulla et al 2019 scripts without an issue
    source('R/Prepare_Data_from_Bulla_et_al_2019.R')# generates 18 warnings, same way as Kubelke et al's  2018 script
    b = data.table(b)
    b[,hem01 := ifelse(hemisphere == 'Northern', 1,0)]

# check whether Beintema transformation of Kubelka follows from their methods (i.e. vs our recalculated ones) - it mostly does
    u=b[b$DPRtrans=='YES',]
    #u=u[-which(is.na(u$obs_time)| is.na(u$other_failed)),]
    r = u$obs_time
    u$expMB = (r*u$Incubation_days*(u$other_failed +u$predated)/2)+(r*u$Incubation_days*(u$hatched+u$infertile))
    u[Exposure_days!=expMB,.(Exposure_days, expMB)]
    u$DPR_MB = u$predated/u$expMB
    u[DPR_orig!=DPR_MB,.(DPR_orig, DPR_MB)]
    summary(u[DPR_orig!=DPR_MB,DPR_orig-DPR_MB])

# recalculate DPR and add DPR based on 50% assumption for transformed data
    b[, exposure_recal := Exposure_days]
    b[DPRtrans=='YES', exposure_recal := (obs_time*Incubation_days*(other_failed +predated)/2)+(obs_time*Incubation_days*(hatched+infertile))]
    b[, DPR_orig_recalc:= predated/exposure_recal]
    b[is.na(DPR_orig_recalc), DPR_orig_recalc := DPR_orig]

    r = 0.5
    b[, exposure_trans50 := exposure_recal]
    b[, DPRtrans50 := DPR_orig_recalc]
    b[DPRtrans=='YES', exposure_trans50 :=(r*Incubation_days*(other_failed +predated)/2)+(r*Incubation_days*(hatched+infertile))]
    b[DPRtrans=='YES',DPRtrans50:= predated/exposure_trans50]

# explore the differences
    summary(b$DPRtrans50-b$DPR_orig_recalc)

    ggplot(b[DPRtrans=='YES'], aes(x = log(DPR_orig_recalc+0.01), y = log(DPRtrans50+0.01))) + geom_point()
    ggplot(b[DPRtrans=='YES'], aes(x = (DPR_orig_recalc), y = (DPRtrans50))) + geom_point()

    ggplot(b[DPRtrans=='YES'], aes(x = log(DPR_orig_recalc+0.01) - log(DPRtrans50+0.01), y = Latitude)) + geom_point()
    ggplot(b[DPRtrans=='YES'], aes(x = DPR_orig_recalc - DPRtrans50, y = Latitude)) + geom_point()
    #ggplot(b, aes(x = log(DPR_orig_recalc+0.01), y = log(DPRtrans50+0.01))) + geom_point()

# Table - only predictors scaled -  linear and mixed models with original and 50% transformation
  # linear  
    m0 = lm(log(DPR_orig+0.01) ~ log(N_nests) + hemisphere + scale(mean_year)+scale(abs(Latitude)),  data = b)
    summary(glht(m0))
    m0_50 = lm(log(DPRtrans50+0.01) ~ log(N_nests) + hemisphere + scale(mean_year)+scale(abs(Latitude)), data = b)
    summary(glht(m0_50))

    m3 = lm(log(DPR_orig+0.01) ~ log(N_nests) + hemisphere*scale(mean_year)*scale(abs(Latitude)),  data = b)
    summary(glht(m3))
    m3_50 = lm(log(DPRtrans50+0.01) ~ log(N_nests) + hemisphere*scale(mean_year)*scale(abs(Latitude)), data = b)
    summary(glht(m3_50))

     m2 = lm(log(DPR_orig+0.01) ~ log(N_nests) + hemisphere + scale(mean_year)*scale(abs(Latitude)),  data = b)
    summary(glht(m2))
    m2_50 = lm(log(DPRtrans50+0.01) ~ log(N_nests) + hemisphere + scale(mean_year)*scale(abs(Latitude)), data = b)
    summary(glht(m2_50))
  # mixed    
    mm0 = lmer(log(DPR_orig+0.01) ~ log(N_nests) + hemisphere + scale(mean_year)+scale(abs(Latitude)) + (1|site), data = b)
    summary(glht(mm0))

    mm0_50 = lmer(log(DPRtrans50+0.01) ~ log(N_nests) + hemisphere + scale(mean_year)+scale(abs(Latitude)) + (1|site), data = b)
    summary(glht(mm0_50))

    mm3 = lmer(log(DPR_orig+0.01) ~ log(N_nests) + hemisphere*scale(mean_year)*scale(abs(Latitude)) + (1|site), data = b)
    summary(glht(mm3))

    mm3_50 = lmer(log(DPRtrans50+0.01) ~ log(N_nests) + hemisphere*scale(mean_year)*scale(abs(Latitude)) + (1|site), data = b)
    summary(glht(mm3_50))

    mm2 = lmer(log(DPR_orig+0.01) ~ log(N_nests) + hemisphere + scale(mean_year)*scale(abs(Latitude)) + (1|site), data = b)
    summary(glht(mm2))

    mm2_50 = lmer(log(DPRtrans50+0.01) ~ log(N_nests) + hemisphere + scale(mean_year)*scale(abs(Latitude)) + (1|site), data = b)
    summary(glht(mm2_50))
  # export model results - Table
    o1 = m_out_lm(name = "linear Kubelka", model = m3, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
    o2 = m_out_lm(name = "linear 50%", model = m3_50, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
    o3 = m_out(name = "mixed Kubelka", model = mm3, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
    o4 = m_out(name = "mixed 50%", model = mm3_50, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
    
    o1$p_value = NA
    o1$p_value[1:9] = round(summary(glht(m3))$test$pvalues,2)
    o2$p_value = NA
    o2$p_value[1:9] = round(summary(glht(m3_50))$test$pvalues,2)
    o3$p_value = NA
    o3$p_value[1:9] = round(summary(glht(mm3))$test$pvalues,2)
    o4$p_value = NA
    o4$p_value[1:9] = round(summary(glht(mm3_50))$test$pvalues,2)

   write_xlsx(rbind(o1,o2,o3,o4), 'Outputs/Table_Note6.xlsx')

# plot - all scaled
  # prepare model predictions
    # linear  
      m3 = lm(scale(log(DPR_orig+0.01)) ~ scale(log(N_nests)) + scale(hem01)*scale(mean_year)*scale(abs(Latitude)),  data = b)
      summary(glht(m3))
      m3_50 = lm(scale(log(DPR_orig+0.01)) ~ scale(log(N_nests)) + scale(hem01)*scale(mean_year)*scale(abs(Latitude)), data = b)
      summary(glht(m3_50))
    # mixed    
      mm3 = lmer(scale(log(DPR_orig+0.01)) ~ scale(log(N_nests)) + scale(hem01)*scale(mean_year)*scale(abs(Latitude)) + (1|site), data = b)
      summary(glht(mm3))

      mm3_50 = lmer(scale(log(DPR_orig+0.01)) ~ scale(log(N_nests)) + scale(hem01)*scale(mean_year)*scale(abs(Latitude)) + (1|site), data = b)
      summary(glht(mm3_50))

    # export results 
      o1 = m_out_lm(name = "linear Kubelka", model = m3, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
      o2 = m_out_lm(name = "linear 50%", model = m3_50, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
      o3 = m_out(name = "mixed Kubelka", model = mm3, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
      o4 = m_out(name = "mixed 50%", model = mm3_50, round_ = 3, nsim = 5000, aic = FALSE)[1:6]
  # plot
    o =data.table(rbind(o1,o2,o3,o4))
    o = o[type == 'fixed' & !effect%in%c('(Intercept)')]
    o[, estimate_r := as.numeric(estimate_r)]
    o[, lwr_r := as.numeric(lwr_r)]
    o[, upr_r := as.numeric(upr_r)]

    o[, effect2 := factor(effect, levels = c(
      "scale(log(N_nests))",
      "scale(hem01)", 
      "scale(mean_year)",
      "scale(abs(Latitude))",  
      "scale(hem01):scale(mean_year)", 
      "scale(hem01):scale(abs(Latitude))",
      "scale(mean_year):scale(abs(Latitude))",
      "scale(hem01):scale(mean_year):scale(abs(Latitude))"))]

    o[, effect3 := effect2]
    o[ effect3 %in% c("scale(log(N_nests))"), effect3 := "N nests"]
    o[ effect3 %in% c("scale(hem01)"), effect3 := "Hemisphere"]
    o[ effect3 %in% c("scale(mean_year)"), effect3 := "Mean year"]
    o[ effect3 %in% c("scale(abs(Latitude))"), effect3 := "Latitude(absolute)"]
    o[ effect3 %in% c("scale(hem01):scale(mean_year)"), effect3 := "Hemisphere x Mean year"]
    o[ effect3 %in% c("scale(hem01):scale(abs(Latitude))"), effect3 := "Hemisphere x Latitude"]
    o[ effect3 %in% c("scale(mean_year):scale(abs(Latitude))"), effect3 := "Mean year x Latitude"]
    o[ effect3 %in% c("scale(hem01):scale(mean_year):scale(abs(Latitude))"), effect3 := "Hemisphere x Mean year x Latitude"]
  
    o[, model2:= model]
    o[ model2 %in% c("linear Kubelka"), model2 := "linear | as in Kubelka et al 2018"]
    o[ model2 %in% c("linear 50%"), model2 := "linear | 50% of incubation period"]
    o[ model2 %in% c("mixed Kubelka"), model2 := "mixed | as in Kubelka et al 2018"]
    o[ model2 %in% c("mixed 50%"), model2 := "mixed | 50% of incubation period"]
    o[,model2:= factor(model2, levels = c(
                "mixed | 50% of incubation period",
                "mixed | as in Kubelka et al 2018",
                "linear | 50% of incubation period",
                "linear | as in Kubelka et al 2018" 
                ) )]

    g=
    ggplot(o, aes( x = estimate_r, y = effect3, col = model2)) +
      geom_vline(xintercept = 0, color = 'grey',linetype = "dotted") +
      geom_errorbarh(aes(xmin=lwr_r, xmax=upr_r), height=0, position=ggstance::position_dodgev(.6)) +  
      geom_point(position=ggstance::position_dodgev(.6))+
      ylab('') + xlab('Standardized effect size') +
      scale_colour_brewer(type = 'qual', palette = 'Paired', name = "Model | DPR estimation for unknown",guide = guide_legend(reverse = TRUE))+#labels = unique(o$model2)[4:1]  , breaks = as.character(unique(o$model2)[4:1])) +# labels = legend_label[12:1], breaks = as.character(nn$set2[12:1]) ) 
      scale_x_continuous(limits = c(-0.27, 0.66), expand = c(0, 0), breaks = seq(-0.2,0.6, by = 0.2), labels = seq(-0.2,0.6, by = 0.2)) +
      theme_bw()+
      theme(      legend.position ="right",
                  legend.title=element_text(size=7), 
                  legend.text=element_text(size=6),
                  #legend.spacing.y = unit(0.1, 'cm'), 
                  legend.key.height= unit(0.5,"line"),
                  plot.margin = margin(b = 0.5, l = 0.5, t = 0.5, r =0.5, unit =  "pt"),
                  panel.grid = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(colour = ax_lines, size = 0.25),
                  axis.line.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.ticks.x= element_line( colour = ax_lines, size = 0.25),
                  #axis.text.x = element_text()
                  axis.ticks.length = unit(1, "pt"),
                  axis.text.x=element_text(, size = 6),
                  axis.text.y=element_text(colour="black", size = 7),
                  axis.title=element_text(size=7)
                  )

   ggsave(file="Outputs/Figure_Note6.png",g, dpi = 600, width = 18, height = 4.5, units = "cm") 

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
  #  [1] readxl_1.3.1       ape_5.4-1          writexl_1.3.1      viridis_0.5.1     
  #  [5] viridisLite_0.3.0  RColorBrewer_1.1-2 performance_0.4.8  multcomp_1.4-13   
  #  [9] TH.data_1.0-10     survival_3.1-12    mvtnorm_1.1-1      ggplot2_3.3.2     
  # [13] data.table_1.13.0  arm_1.11-2         lme4_1.1-23        Matrix_1.2-18     
  # [17] MASS_7.3-51.6     
  # 
  # loaded via a namespace (and not attached):
  #  [1] Rcpp_1.0.5          here_0.1            lattice_0.20-41     png_0.1-7          
  #  [5] zoo_1.8-8           rprojroot_1.3-2     digest_0.6.25       cellranger_1.1.0   
  #  [9] R6_2.4.1            backports_1.1.8     acepack_1.4.1       coda_0.19-3        
  # [13] pillar_1.4.6        rlang_0.4.7         rstudioapi_0.11     minqa_1.2.4        
  # [17] nloptr_1.2.2.2      rpart_4.1-15        checkmate_2.0.0     splines_4.0.2      
  # [21] statmod_1.4.34      stringr_1.4.0       foreign_0.8-80      htmlwidgets_1.5.1  
  # [25] munsell_0.5.0       compiler_4.0.2      xfun_0.16           pkgconfig_2.0.3    
  # [29] base64enc_0.1-3     htmltools_0.5.0     nnet_7.3-14         insight_0.9.0      
  # [33] tidyselect_1.1.0    tibble_3.0.3        gridExtra_2.3       htmlTable_2.0.1    
  # [37] Hmisc_4.4-0         codetools_0.2-16    crayon_1.3.4        dplyr_1.0.1        
  # [41] withr_2.2.0         grid_4.0.2          nlme_3.1-148        gtable_0.3.0       
  # [45] lifecycle_0.2.0     magrittr_1.5        bayestestR_0.7.2    scales_1.1.1       
  # [49] stringi_1.5.3       latticeExtra_0.6-29 ellipsis_0.3.1      generics_0.0.2     
  # [53] vctrs_0.3.2         boot_1.3-25         sandwich_2.5-1      Formula_1.2-3      
  # [57] tools_4.0.2         glue_1.4.2          purrr_0.3.4         jpeg_0.1-8.1       
  # [61] parallel_4.0.2      abind_1.4-5         colorspace_1.4-1    cluster_2.1.0      
  # [65] knitr_1.29 