require(ggplot2)

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
