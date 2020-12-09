require(ggplot2)
require(data.table)

data <- read.csv( "Data/DATApopulations_research intensity.csv", h = T, sep=";")
data$Belt = factor(data$Belt, levels = c("Arctic","North temperate", "North tropics","South tropics","South temperate"))

ggplot(data[!is.na(data$research_intensity),], aes(y = mean_year, x = research_intensity)) + 
  geom_violin(trim = TRUE, fill = '#A4A4A4', col = NA) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +  
  geom_boxplot(width=0.1, col = 'red', outlier.size = 1, outlier.shape = NA,fill = NA)+
  
  facet_wrap(~Belt, ncol = 3, dir = "v") +
  coord_flip() +

  scale_x_discrete("Research intensity", labels = c("small", "medium", "hight"))+
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

# update 
  d = data.table(data)
  d[mean_year<1960, period := '1943-1959']
  d[is.na(period) & mean_year<1970, period := '1960-1969']
  d[is.na(period) & mean_year<1980, period := '1970-1979']
  d[is.na(period) & mean_year<1990, period := '1980-1989']
  d[is.na(period) & mean_year<2000, period := '1990-1999']
  d[is.na(period) & mean_year<2010, period := '2000-2009']
  d[is.na(period) & mean_year<2020, period := '2010-2016']

  d[mean_year<1960, period2 := '1943-1959']
  d[is.na(period2) & mean_year<1980, period2 := '1960-1979']
  d[is.na(period2) & mean_year<2000, period2 := '1980-1999']
  d[is.na(period2) & mean_year<2020, period2 := '2000-2016']

  # count like Kubelka
    ggplot(d[!is.na(research_intensity)], aes(x = period2, fill = research_intensity)) + geom_bar() + 
    facet_wrap(~Belt, ncol =1) + xlab("Period") + ylab("Populations [count]") +
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
  # count detailed
    ggplot(d[!is.na(research_intensity)], aes(x = period, fill = research_intensity)) + geom_bar() + 
    facet_wrap(~Belt, ncol =1) + xlab("Period") + ylab("Populations [count]") +
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

      ggsave(file="Outputs/temp_change_res-int_barCount_detailed.png", dpi = 600, width = 8, height = 11, units = "cm")
  # percentage like Kubelka
    dd = d[, .(populations = .N), by = list(period2, Belt, research_intensity)]
    dd[, n := sum(populations), by = list(period2, Belt)]
    ggplot(dd[!is.na(research_intensity)], aes(x = period2, y= 100*populations/n, fill = research_intensity)) + 
    geom_bar(position = "fill",stat = "identity") +
    scale_y_continuous(labels = scales::percent_format()) +
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

  # percentage detailed  
    dd = d[, .(populations = .N), by = list(period, Belt, research_intensity)]
    dd[, n := sum(populations), by = list(period, Belt)]
    ggplot(dd[!is.na(research_intensity)], aes(x = period, y= 100*populations/n, fill = research_intensity)) + 
      geom_bar(position = "fill",stat = "identity") +
      scale_y_continuous(labels = scales::percent_format()) +
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

    ggsave(file="Outputs/temp_change_res-int_barPerc_detailed.png", dpi = 600, width = 8, height = 11, units = "cm")
  