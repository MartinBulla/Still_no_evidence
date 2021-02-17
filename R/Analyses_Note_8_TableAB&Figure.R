# ===========================================================================
# Supporting information for "Bulla, Valcu & Kempenaers. (2021) "Still no 
# evidence for disruption of global patterns of nest predation in shorebirds" 
# Contributor: Martin Bulla
# 📍 script runs relative to the project's root directory and uses Kubelka et 
# al's and Bulla et al.'s 2019 data to generate Note's 9 Table AB and Figure
# ===========================================================================

# TOOLS
  rm( list = ls() ) 
  PNG = FALSE # print figures in PNG or not
  source('R/Constants_Functions.R')
  require(plyr); require(magrittr); require(ggnewscale)

# data preparation  
    d = data.table(read.csv("Data/Rebuttal_data_with_DPR_by_yr.csv", h = T, sep=",",stringsAsFactors = FALSE))
    #dd = d[year!='all' & duplicated(paste(source_id,species, locality, year)), .(source_id, species, locality,year)] 
    #print(d[paste(source_id,species, locality, year)%in%dd[,paste(source_id,species, locality, year)],.(who, source_id,species, locality, Latitude, year, N.nests)], nrows =200)

    d[,pk := 1:nrow(d)]
    d[,Belt := as.factor(Belt)]
    d[,abs_lat := abs(Latitude)]
    d[,hemisphere :=as.factor(ifelse(Latitude > 0, "Northern", "Southern"))]
    d[,site := paste(Latitude, Longitude)]
    d[,year_:=as.numeric(ifelse(year=='all', ((start_year+end_year)/2), d$year))]
    d = d[-which(is.na(year)|is.na(DPR))]
    d = d[which(d$year!='all')]
    summary(factor(d$DPR_assumption))
    d[is.na(DPR_assumption), DPR_assumption:='our_est']
    d[, N_nests_loc_year := sum(N.nests), by =list(species, locality, year)] # localities are larger or there are estimates for various habitats, so there are sometimes more estimates per locality in a given year
    
    b = d[N.nests>11] # uses more than 11 nests for fine scaled localities and the output is same as for d = d[d$N_nests_loc_year>11,] (except for Vanellus_vanellus source_id 118 where additional estimates based on 6,8, and 9 would be included) 

    summary(factor(b$DPR_assumption))
    b[,fail_is_depr := ifelse(DPR_assumption == 'our_est', 'No','Yes')]

    b$year = as.numeric(b$year)
    b$n=1
    b$pop = paste(b$source_id,b$species,b$locality)
    b$direct = ifelse(b$DPRtrans == "YES", "No", "Yes")
    b= b[!(locality =="Barrow, Alaska" & year%in%c( 2003,2004))]    
    x = ddply(b,.(pop), summarise, n =length(unique(year)))
    summary(factor(x$n))
    nrow(x[x$n>5,])
    #x2 = ddply(b,.(pop), summarise, n =sum(n))
    #nrow(x2[x2$n>5,])
    #d[!pk%in%b$pk[b$species=='Vanellus_vanellus'] & species=='Vanellus_vanellus']

    o = data.table(b[b$pop%in%x$pop[x$n>5],])
    od = o[o$fail_is_depr == "No",]
    
    o[,  slope := lm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_dir := ifelse(slope >0, '+', '-')]
    o[,  slope_rlm_z := rlm(scale(DPR) ~ scale(year_))  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_rlm := rlm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_dir_rlm := ifelse(slope_rlm >0, '+', '-')]
    o[,  slope_rlm_p := summary(glht(rlm(scale(DPR) ~ scale(year_))))$test$pvalues[2], by = pop]
    o[, labeler := paste(round(Latitude), "|", species)]
    o[slope_rlm_p<0.1, labeler :=paste(labeler,'*')]
    o[slope_rlm_p<0.01, labeler :=paste0(labeler,'*')]
    o[slope_rlm_p<0.001,labeler := paste0(labeler,'*')]
    o[, labeler := factor(labeler)]
    o[, labeler := factor(labeler, levels=rev(levels(labeler)))]

    oo =o[!duplicated(labeler)]
    oo[slope_rlm_p<0.25, round(slope_rlm_p,4)]
    oo[order(Latitude) & slope_rlm_p<0.25,.(Latitude,species,locality,slope_rlm_p)]
    summary(factor(oo$slope_dir))
    summary(factor(oo$slope_dir[oo$fail_is_depr =='No']))
    xtabs(~oo$slope_dir+oo$fail_is_depr+oo$direct) 
    summary(oo$slope_rlm_z)
    t.test(oo$slope_rlm[oo$Belt == 'Arctic'], oo$slope_rlm[oo$Belt != 'Arctic'])  
  
# TABLE 
    nrow(o); length(unique(o$year_)); length(unique(o$site)); length(unique(o$species))
    nrow(od); length(unique(od$year_)); length(unique(od$site)); length(unique(od$species))

      xx = ddply(o,.(pop, locality), summarise, n =length(unique(year)))
      xy = ddply(xx,.(locality), summarise, n =sum(n))
    
    # random slope within site
    m1sl = lmer(log(DPR+0.01)~scale(log(N.nests))+scale(Latitude) + scale(year_) +(1|species) + (scale(year_)|locality),o)
    summary(m1sl)
    summary(glht(m1sl))

    # random slope within species - same results
    #m1s = lmer(log(DPR+0.01)~scale(log(N.nests))+scale(Latitude) + scale(year_) +(scale(year_)|species) + (1|locality),o)
    #summary(m1s)
    #summary(glht(m1sl))

    # without random slope - same results
    #m1 = lmer(log(DPR+0.01)~scale(log(N.nests))+scale(Latitude) + scale(year_) +(1|species) + (1|locality),o)
    #summary(m1)
    #summary(glht(m1))

    # controlling for failure=predation cases - same results
    #mf = lmer(log(DPR+0.01)~scale(log(N.nests))+fail_is_depr+scale(Latitude) + scale(year_) +(1|species) + (scale(year_)|locality),o)
    #mf = lmer(log(DPR+0.01)~scale(log(N.nests))+fail_is_depr+scale(Latitude) + fail_is_depr*scale(year_) +(1|species) + (scale(year_)|locality),o)
    #summary(mf)
    #summary(glht(mf))

    m2sl = lmer(log(DPR+0.01)~scale(log(N.nests))+scale(Latitude) + scale(year_) +(1|species) + (scale(year_)|locality),od)
    summary(m2sl)
    summary(glht(m2sl))

    #m2s = lmer(log(DPR+0.01)~scale(log(N.nests))+scale(Latitude) + scale(year_) +(scale(year_)|species) + (1|locality),od)
    #summary(m2s)
    #summary(glht(m2s))

    #m2 = lmer(log(DPR+0.01)~scale(log(N.nests))+scale(Latitude) + scale(year_) +(1|species) + (1|locality),od)
    #summary(m2)
    #summary(glht(m2))
    
    o1 = m_out(name = "all", model = m1sl, round_ = 3, nsim = 5000, aic = FALSE)  
    o2 = m_out(name = "pred_only", model = m2sl, round_ = 3, nsim = 5000, aic = FALSE)  
    
    o1$p_value = NA
    o1$p_value[1:4] = round(summary(glht(m1sl))$test$pvalues,2)
    o2$p_value = NA
    o2$p_value[1:4] = round(summary(glht(m2sl))$test$pvalues,2)
    

    sname = 'Table_Note8'
    write_xlsx(rbind(o1,o2), paste0("Outputs/",sname,'.xlsx'))
    
# PLOT
    direct = "#FCB42C" #"#5eab2b"#c3dbe5" #9CC3D5""#FCB42C"   #ffa500
    trans = "#535F7C" #"#7c9caa" #c3dbe5" #9CC3D5##7ac5cd"# "#535F7C"  #0063B2FF" #"#619CFF" "#9CC3D5FF"
    fam_lines_pos =  "grey80"##FCB42C" 
    fam_lines_neg = "grey30"##535F7C"
    all_failed = 17
    pred_only = 16
    ggplot(o, aes(x = year_, y = DPR)) + 
          geom_smooth(aes(color = slope_dir_rlm, fill = slope_dir_rlm), method = 'rlm',  size = 0.5) +
          scale_color_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
          scale_fill_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
          geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +

          new_scale_color() +   
          geom_point(aes(col = direct, shape = fail_is_depr), size = 1.1) + 
          scale_color_manual(name = "Directly\ncalculated\nestimates:",
                             values = c(direct, trans),
                             breaks = c("Yes", "No"))+
          scale_shape_manual(name = "All failed\nassumed as\npredated:",
                             values = c(all_failed, pred_only),
                             breaks = c("Yes", "No"))+

          facet_wrap(~labeler, ncol = 3, scale = "free_y", as.table = TRUE)+ 
          ylab("Daily predation rate") + xlab("Year") +
          #labs(col = "Directly calculated:")+
          
          theme_light()+
                          theme(  axis.line=element_line(colour="grey70", size=0.25),
                                  #panel.border=element_rect(colour="white"),
                                  panel.border=element_rect(colour="grey70", size=0.25),
                                  panel.grid = element_blank(),
                                  
                                  axis.title=element_text(size=7, colour="grey30"),
                                  axis.title.y = element_text(vjust=1),
                                  axis.title.x = element_text(vjust=1),
                                  axis.text=element_text(size=6),# margin=units(0.5,"mm")),
                                  axis.ticks.length=unit(0.5,"mm"),
                                  #axis.ticks.margin,
                                  
                                  strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
                                  strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
                                  #strip.background = element_blank(), 
                                  #strip.text = element_blank(),
                                  panel.spacing = unit(0, "mm"),
                                  #legend.position="none"
                                  
                                  legend.key=element_rect(fill="grey99", colour="white"),
                                  legend.text=element_text(size=6, colour="grey30"),
                                  legend.title=element_text(size=6, colour="grey30"),
                                  legend.key.height= unit(0.5,"line"),
                                  legend.key.width = unit(0.25, "cm"),
                                  legend.margin = margin(0,0,0,0, unit="cm"),
                                  legend.box.margin = margin(l = -5),
                                  legend.background = element_blank()
                                  #legend.justification = c(-1,0),
                                  #legend.spacing.y = unit(.5, "cm")
                                  )

    ggsave(file="Outputs/Figure_Note8.png", dpi = 600, width = 14, height = 17, unit = "cm")

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
  #  [1] ggnewscale_0.4.3   magrittr_1.5       plyr_1.8.6         writexl_1.3.1     
  #  [5] viridis_0.5.1      viridisLite_0.3.0  RColorBrewer_1.1-2 performance_0.4.8 
  #  [9] multcomp_1.4-13    TH.data_1.0-10     survival_3.1-12    mvtnorm_1.1-1     
  # [13] here_0.1           ggplot2_3.3.2      data.table_1.13.0  arm_1.11-2        
  # [17] lme4_1.1-23        Matrix_1.2-18      MASS_7.3-51.6     
  # 
  # loaded via a namespace (and not attached):
  #  [1] Rcpp_1.0.5          lattice_0.20-41     png_0.1-7           zoo_1.8-8          
  #  [5] rprojroot_1.3-2     digest_0.6.25       R6_2.4.1            backports_1.1.8    
  #  [9] acepack_1.4.1       coda_0.19-3         pillar_1.4.6        rlang_0.4.7        
  # [13] rstudioapi_0.11     minqa_1.2.4         nloptr_1.2.2.2      rpart_4.1-15       
  # [17] checkmate_2.0.0     labeling_0.3        splines_4.0.2       statmod_1.4.34     
  # [21] stringr_1.4.0       foreign_0.8-80      htmlwidgets_1.5.1   munsell_0.5.0      
  # [25] compiler_4.0.2      xfun_0.16           pkgconfig_2.0.3     base64enc_0.1-3    
  # [29] mgcv_1.8-31         htmltools_0.5.0     insight_0.9.0       nnet_7.3-14        
  # [33] tidyselect_1.1.0    tibble_3.0.3        gridExtra_2.3       htmlTable_2.0.1    
  # [37] Hmisc_4.4-0         codetools_0.2-16    crayon_1.3.4        dplyr_1.0.1        
  # [41] withr_2.2.0         grid_4.0.2          nlme_3.1-148        gtable_0.3.0       
  # [45] lifecycle_0.2.0     bayestestR_0.7.2    scales_1.1.1        stringi_1.5.3      
  # [49] farver_2.0.3        latticeExtra_0.6-29 ellipsis_0.3.1      generics_0.0.2     
  # [53] vctrs_0.3.2         boot_1.3-25         sandwich_2.5-1      Formula_1.2-3      
  # [57] tools_4.0.2         glue_1.4.2          purrr_0.3.4         jpeg_0.1-8.1       
  # [61] abind_1.4-5         colorspace_1.4-1    cluster_2.1.0       knitr_1.29     