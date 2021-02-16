# ===========================================================================
# Supporting information for "Bulla, Valcu & Kempenaers. (2021) "Still no 
# evidence for disruption of global patterns of nest predation in shorebirds" 
# Contributor: Martin Bulla
# 📍 script runs relative to the project's root directory and uses Kubelka et 
# al's data to generate Figure R3 and its variants
# (3) Do predation rates increase over time? 
# ===========================================================================

# load tool and data
  source('R/Constants_Functions.R')
  data <- read.csv( "Data/DATApopulations_research intensity.csv", h = T, sep=";")
  data$Belt = factor(data$Belt, levels = c("Arctic","North temperate", "North tropics","South tropics","South temperate"))
  d = data.table(data)
  d[research_intensity == "3_high",r_int:="high"]
  d[research_intensity == "2_medium",r_int:="medium"]
  d[research_intensity == "1_small",r_int:="small"]

  # create decades
  d[mean_year<1960, period := '1943-1959']
  d[is.na(period) & mean_year<1970, period := '1960-1969']
  d[is.na(period) & mean_year<1980, period := '1970-1979']
  d[is.na(period) & mean_year<1990, period := '1980-1989']
  d[is.na(period) & mean_year<2000, period := '1990-1999']
  d[is.na(period) & mean_year<2010, period := '2000-2009']
  d[is.na(period) & mean_year<2020, period := '2010-2016']

  # create larger periods
  d[mean_year<1960, period2 := '1943-1959']
  d[is.na(period2) & mean_year<1980, period2 := '1960-1979']
  d[is.na(period2) & mean_year<2000, period2 := '1980-1999']
  d[is.na(period2) & mean_year<2020, period2 := '2000-2016']

# plot as counts and detailed
    ggplot(d[!is.na(research_intensity)], aes(x = period, fill = r_int)) + geom_bar() + 
    facet_wrap(~Belt, ncol =1) + xlab("Period") + ylab("Populations [count]") +
    scale_fill_viridis(discrete=TRUE, name = "Research intensity") +
    theme_light()+ theme(  
          axis.line=element_line(colour="grey70", size=0.25),
          #panel.border=element_rect(colour="white"),
          panel.border=element_rect(colour="grey70", size=0.25),
          panel.grid = element_blank(),
          
          axis.title=element_text(size=7, colour="grey30"),
          axis.title.y = element_text(vjust=1),
          axis.title.x = element_text(vjust=1),
          axis.text=element_text(size=6),# margin=units(0.5,"mm")),
          axis.text.x=element_text(angle = 45, hjust=1),
          axis.ticks.length=unit(0.5,"mm"),
          #axis.ticks.margin,
          
          strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
          strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
          #strip.background = element_blank(), 
          #strip.text = element_blank(),
          legend.title = element_text(size=6, margin = margin(l = -4, unit = "pt")),
          legend.text=element_text(size=6, margin = margin(l = -2, unit = "pt")),
          legend.key.height= unit(0.5,"line"),
          legend.key.width = unit(0.2, "cm"),
          panel.spacing = unit(0, "mm")
          #legend.position="none"
          )

      ggsave(file="Outputs/Figure_R3left.png", dpi = 600, width = 8, height = 11, units = "cm")
# plot as percentages and detailed  
    dd = d[, .(populations = .N), by = list(period, Belt, r_int)]
    dd[, n := sum(populations), by = list(period, Belt)]
    ggplot(dd[!is.na(r_int)], aes(x = period, y= 100*populations/n, fill = r_int)) + 
      geom_bar(position = "fill",stat = "identity") +
      scale_y_continuous(labels = scales::percent_format()) +
      scale_fill_viridis(discrete=TRUE, name = "Research intensity") +
      facet_wrap(~Belt, ncol =1) + xlab("Period") + ylab("Populations") +
      theme_light()+ theme(  
        axis.line=element_line(colour="grey70", size=0.25),
        #panel.border=element_rect(colour="white"),
        panel.border=element_rect(colour="grey70", size=0.25),
        panel.grid = element_blank(),
        
        axis.title=element_text(size=7, colour="grey30"),
        axis.title.y = element_text(vjust=1),
        axis.title.x = element_text(vjust=1),
        axis.text=element_text(size=6),# margin=units(0.5,"mm")),
        axis.text.x=element_text(angle = 45, hjust=1),
        axis.ticks.length=unit(0.5,"mm"),
        #axis.ticks.margin,
        
        strip.text.x =element_text(size = 6, color="grey30",  margin=margin(1,1,1,1,"mm")), #grey50
        strip.background=element_rect(fill="grey99",colour="grey70", size=0.25),
        #strip.background = element_blank(), 
        #strip.text = element_blank(),
        legend.title = element_text(size=6, margin = margin(l = -4, unit = "pt")),
        legend.text=element_text(size=6, margin = margin(l = -2, unit = "pt")),
        legend.key.height= unit(0.5,"line"),
        legend.key.width = unit(0.2, "cm"),
        panel.spacing = unit(0, "mm")
        #legend.position="none"
        )

    ggsave(file="Outputs/Figure_R3_right.png", dpi = 600, width = 8, height = 11, units = "cm")
  
# not used - count like Kubelka et al period division, but per Belt
    ggplot(d[!is.na(r_int)], aes(x = period2, fill = r_int)) + geom_bar() + 
    facet_wrap(~Belt, ncol =1) + xlab("Period") + ylab("Populations [count]") +
    scale_fill_viridis(discrete=TRUE, name = "Research intensity") +
    theme_light()+ theme(  
          axis.line=element_line(colour="grey70", size=0.25),
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
          legend.title = element_text(size=6, margin = margin(l = -4, unit = "pt")),
          legend.text=element_text(size=6, margin = margin(l = -2, unit = "pt")),
          legend.key.height= unit(0.5,"line"),
          legend.key.width = unit(0.2, "cm"),
          panel.spacing = unit(0, "mm")
          #legend.position="none"
          )

      ggsave(file="Outputs/temp_change_res-int_barCount.png", dpi = 600, width = 8, height = 11, units = "cm")    
# not used - percentage like Kubelka et al period division, but per Belt
    dd = d[, .(populations = .N), by = list(period2, Belt, r_int)]
    dd[, n := sum(populations), by = list(period2, Belt)]

    ggplot(dd[!is.na(r_int)], aes(x = period2, y= 100*populations/n, fill = r_int)) + 
    geom_bar(position = "fill",stat = "identity") +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_viridis(discrete=TRUE, name = "Research intensity") +
      facet_wrap(~Belt, ncol =1) + xlab("Period") + ylab("Populations") +
      theme_light()+ theme(  
        axis.line=element_line(colour="grey70", size=0.25),
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
        legend.title = element_text(size=6, margin = margin(l = -4, unit = "pt")),
        legend.text=element_text(size=6, margin = margin(l = -2, unit = "pt")),
        legend.key.height= unit(0.5,"line"),
        legend.key.width = unit(0.2, "cm"),
        panel.spacing = unit(0, "mm")
        #legend.position="none"
        )

      ggsave(file="Outputs/temp_change_res-int_barPerc.png", dpi = 600, width = 8, height = 11, units = "cm")
# not used - initial figure from the text sent to Kubelka et al.
  ggplot(data[!is.na(data$research_intensity),], aes(y = mean_year, x = research_intensity)) + 
    geom_violin(trim = TRUE, fill = '#A4A4A4', col = NA) + 
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +  
    geom_boxplot(width=0.1, col = 'red', outlier.size = 1, outlier.shape = NA,fill = NA)+
    
    facet_wrap(~Belt, ncol = 3, dir = "v") +
    coord_flip() +

    scale_x_discrete("Research intensity", labels = c("small", "medium", "high"))+
    scale_y_continuous("Mean year", limits = c(1935, 2020), breaks = c(1940, 1960, 1980,2000,2020))+

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
                                  panel.spacing = unit(0, "mm")
                                  #legend.position="none"
                                  )

  ggsave(file="Outputs/temp_change_res-int_3col.png", dpi = 600, width = 12, height = 8, units = "cm")

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
  #  [1] writexl_1.3.1      viridis_0.5.1      viridisLite_0.3.0  RColorBrewer_1.1-2
  #  [5] performance_0.4.8  multcomp_1.4-13    TH.data_1.0-10     survival_3.1-12   
  #  [9] mvtnorm_1.1-1      ggplot2_3.3.2      data.table_1.13.0  arm_1.11-2        
  # [13] lme4_1.1-23        Matrix_1.2-18      MASS_7.3-51.6     
  # 
  # loaded via a namespace (and not attached):
  #  [1] Rcpp_1.0.5          lattice_0.20-41     png_0.1-7           zoo_1.8-8          
  #  [5] digest_0.6.25       R6_2.4.1            backports_1.1.8     acepack_1.4.1      
  #  [9] coda_0.19-3         pillar_1.4.6        rlang_0.4.7         rstudioapi_0.11    
  # [13] minqa_1.2.4         nloptr_1.2.2.2      rpart_4.1-15        checkmate_2.0.0    
  # [17] labeling_0.3        splines_4.0.2       statmod_1.4.34      stringr_1.4.0      
  # [21] foreign_0.8-80      htmlwidgets_1.5.1   munsell_0.5.0       compiler_4.0.2     
  # [25] xfun_0.16           pkgconfig_2.0.3     base64enc_0.1-3     htmltools_0.5.0    
  # [29] nnet_7.3-14         insight_0.9.0       tidyselect_1.1.0    tibble_3.0.3       
  # [33] gridExtra_2.3       htmlTable_2.0.1     Hmisc_4.4-0         codetools_0.2-16   
  # [37] crayon_1.3.4        dplyr_1.0.1         withr_2.2.0         grid_4.0.2         
  # [41] nlme_3.1-148        gtable_0.3.0        lifecycle_0.2.0     magrittr_1.5       
  # [45] bayestestR_0.7.2    scales_1.1.1        stringi_1.5.3       farver_2.0.3       
  # [49] latticeExtra_0.6-29 ellipsis_0.3.1      generics_0.0.2      vctrs_0.3.2        
  # [53] boot_1.3-25         sandwich_2.5-1      Formula_1.2-3       tools_4.0.2        
  # [57] glue_1.4.2          purrr_0.3.4         jpeg_0.1-8.1        abind_1.4-5        
  # [61] colorspace_1.4-1    cluster_2.1.0       knitr_1.29    