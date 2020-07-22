# TOOLS
  PNG = FALSE # print figures in PNG or not
  rm( list = ls() )	
  sapply(c('AICcmodavg','ape','arm','coxme','data.table','effects', 'ggnewscale','ggplot2','grid', 'lattice','magrittr','mgcv','multcomp','phytools','plyr','RColorBrewer','readxl','writexl'), #XLConnect
      function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE) ))
  source('R/Constants_Functions.R')
# DATA preparation
  d = read.csv("Data/Rebuttal_data_with_DPR_by_yr.csv", h = T, sep=",",stringsAsFactors = FALSE)
  d$pk = 1:nrow(d)
  d$Belt = as.factor(d$Belt)
  d$abs_lat = abs(d$Latitude)
  d$hemisphere =as.factor(ifelse(d$Latitude > 0, "Northern", "Southern"))
  d = d[- which(is.na(d$year)|is.na(d$DPR)),]
  d[is.na(d$DPR),]
  d$year_=as.numeric(ifelse(d$year=='all', ((d$start_year+d$end_year)/2), d$year))
  d = d[d$N.nests>11,]

  b = d[which(d$year!='all'),]
  b$year = as.numeric(b$year)
  b$n=1
  b$pop = paste(b$source_id,b$species,b$locality)
  b$direct = ifelse(b$DPRtrans == "YES", "No", "Yes")
  x = ddply(b,.(pop), summarise, n =sum(n))
  summary(factor(x$n))
  
  o = data.table(b[b$pop%in%x$pop[x$n>5],])
  o[,  slope := lm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
  o[,  slope_dir := ifelse(slope >0, '+', '-')]
  o[,  slope_rlm := rlm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
  o[,  slope_dir_rlm := ifelse(slope_rlm >0, '+', '-')]
  o[, labeler := factor(paste(round(Latitude), "|", species))]
  o[, labeler := factor(labeler, levels=rev(levels(labeler)))]

# TABLE
  bb= b[b$pop%in%x$pop[x$n>1],]
  bb$years_nr = x$n[match(bb$pop,x$pop)]
  nrow(bb)
  length(unique(bb$pop))

  bd = b[b$DPRtrans=='NO',]
  bbd = bb[bb$DPRtrans=='NO',]

  nrow(b); length(unique(b$year_)); length(unique(b$site)); length(unique(b$species))
  nrow(bd); length(unique(bd$year_)); length(unique(bd$locality)); length(unique(bd$species))
  nrow(bb); length(unique(bb$year_)); length(unique(bb$locality)); length(unique(bb$species))
  nrow(bbd); length(unique(bbd$year_)); length(unique(bbd$locality)); length(unique(bbd$species))

  m1 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),b)
  m2 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bd)
  m3 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bb)
  m4 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bbd)
      #summary(glht(m))
      #plot(allEffects(m)) 
  
  o1 = m_out(name = "all", model = m1, round_ = 5, nsim = 5000, aic = FALSE)  
  o2 = m_out(name = "good", model = m2, round_ = 3, nsim = 5000, aic = FALSE)  
  o3 = m_out(name = "all >1yr", model = m3, round_ = 3, nsim = 5000, aic = FALSE)  
  o4 = m_out(name = "good >1yr", model = m4, round_ = 3, nsim = 5000, aic = FALSE)  

  sname = 'Table_Note9'
  write_xlsx(rbind(o1,o2,o3,o4), paste0("Outputs/",sname,'.xlsx'))

# PLOT
  direct = "#FCB42C" #"#5eab2b"#c3dbe5" #9CC3D5""#FCB42C"   #ffa500
  trans = "#535F7C" #"#7c9caa" #c3dbe5" #9CC3D5##7ac5cd"# "#535F7C"  #0063B2FF" #"#619CFF" "#9CC3D5FF"
  fam_lines_pos =  "grey80"##FCB42C" 
  fam_lines_neg = "grey30"##535F7C"
  ggplot(o, aes(x = year_, y = DPR)) + 
        geom_smooth(aes(color = slope_dir_rlm, fill = slope_dir_rlm), method = 'rlm',  size = 0.5) +
        scale_color_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
        scale_fill_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  

        geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +

        new_scale_color() +   
        geom_point(aes(col = direct), pch =16, size = 1.1) + 
        scale_color_manual(name = "Directly\ncalculated\nestimates:",
                           values = c(direct, trans),
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

  ggsave(file="Outputs/rlm_Yearly_pop_points-lines-col.png", dpi = 600, width = 14, height = 17, unit = "cm")

# PLOT alternatives 
  # with RLM lines 
    # lines colored, points colored
      direct = "#FCB42C" #"#5eab2b"#c3dbe5" #9CC3D5""#FCB42C"   #ffa500
      trans = "#535F7C" #"#7c9caa" #c3dbe5" #9CC3D5##7ac5cd"# "#535F7C"  #0063B2FF" #"#619CFF" "#9CC3D5FF"
      fam_lines_pos =  "grey80"##FCB42C" 
      fam_lines_neg = "grey30"##535F7C"
      ggplot(o, aes(x = year_, y = DPR)) + 
            geom_smooth(aes(color = slope_dir_rlm, fill = slope_dir_rlm), method = 'rlm',  size = 0.5) +
            scale_color_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
            scale_fill_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  

            geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +

            new_scale_color() +   
            geom_point(aes(col = direct), pch =16, size = 1.1) + 
            scale_color_manual(name = "Directly\ncalculated\nestimates:",
                               values = c(direct, trans),
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

      ggsave(file="Outputs/rlm_Yearly_pop_points-lines-col.png", dpi = 600, width = 14, height = 17, unit = "cm")
    # lines colored, points colored - 19
      direct = "#FCB42C" #"#5eab2b"#c3dbe5" #9CC3D5""#FCB42C"   #ffa500
      trans = "#535F7C" #"#7c9caa" #c3dbe5" #9CC3D5##7ac5cd"# "#535F7C"  #0063B2FF" #"#619CFF" "#9CC3D5FF"
      fam_lines_pos =  "grey80"##FCB42C" 
      fam_lines_neg = "grey30"##535F7C"
      ggplot(o, aes(x = year_, y = DPR)) + 
            geom_smooth(aes(color = slope_dir_rlm, fill = slope_dir_rlm), method = 'rlm',  size = 0.5) +
            scale_color_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
            scale_fill_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  

            geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +

            new_scale_color() +   
            geom_point(aes(col = direct), alpha = 0.7, size = 1.1) + 
            scale_color_manual(name = "Directly\ncalculated\nestimates:",
                               values = c(direct, trans),
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

      ggsave(file="Outputs/rlm_Yearly_pop_points-lines-col2.png", dpi = 600, width = 14, height = 17, unit = "cm")
    # lines colored, points circles filled/empty
      direct = 19
      trans = 1
      fam_lines_pos =  "#FCB42C" 
      fam_lines_neg = "#535F7C"
      p_col = "#bc507d"#"grey30" #"#bc507d" #7dbc50"#grey50"##5eab2b"
      
      ggplot(o, aes(x = year_, y = DPR)) + 
            geom_smooth(aes(color = slope_dir_rlm, fill = slope_dir_rlm), method = 'rlm',  size = 0.5) +
            scale_color_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
            scale_fill_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  

            geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +

            new_scale_color() +   
            geom_point(aes(shape = direct), size = 1.1, col = p_col) + 
            scale_shape_manual(name = "Directly\ncalculated\nestimates:",
                               values = c(direct, trans),
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

      ggsave(file="Outputs/rlm_Yearly_pop-points-circ_line-col.png", dpi = 600, width = 14, height = 17, unit = "cm")
    # lines colored, points shape
      direct = 19
      trans = 17
      fam_lines_pos =  "#FCB42C"  
      fam_lines_neg = "#535F7C" 
      p_col = "#bc507d" #7dbc50"#grey50"##5eab2b"
      
      ggplot(o, aes(x = year_, y = DPR)) + 
            geom_smooth(aes(color = slope_dir_rlm, fill = slope_dir_rlm), method = 'rlm',  size = 0.5) +
            scale_color_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
            scale_fill_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  

            geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +

            new_scale_color() +   
            geom_point(aes(shape = direct), size = 1.1, col = p_col) + 
            scale_shape_manual(name = "Directly\ncalculated\nestimates:",
                               values = c(direct, trans),
                               breaks = c("Yes", "No"))+

            facet_wrap(~labeler, ncol = 3, scale = "free_y", as.table = TRUE)+ 
            #
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

      ggsave(file="Outputs/rlm_Yearly_pop_points-shape_line-col.png", dpi = 600, width = 14, height = 17, unit = "cm")
  # with LM lines
    # lines colored, points colored
      direct = "#FCB42C" #"#5eab2b"#c3dbe5" #9CC3D5""#FCB42C"   #ffa500
      trans = "#535F7C" #"#7c9caa" #c3dbe5" #9CC3D5##7ac5cd"# "#535F7C"  #0063B2FF" #"#619CFF" "#9CC3D5FF"
      fam_lines_pos =  "grey80"##FCB42C" 
      fam_lines_neg = "grey30"##535F7C"
      ggplot(o, aes(x = year_, y = DPR)) + 
            geom_smooth(aes(color = slope_dir, fill = slope_dir), method = 'lm',  size = 0.5) +
            scale_color_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
            scale_fill_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  

            geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +

            new_scale_color() +   
            geom_point(aes(col = direct), pch =16, size = 1.1) + 
            scale_color_manual(name = "Directly\ncalculated\nestimates:",
                               values = c(direct, trans),
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

      ggsave(file="Outputs/Yearly_pop_points-lines-col.png", dpi = 600, width = 14, height = 17, unit = "cm")
    # lines colored, points colored - 19
      direct = "#FCB42C" #"#5eab2b"#c3dbe5" #9CC3D5""#FCB42C"   #ffa500
      trans = "#535F7C" #"#7c9caa" #c3dbe5" #9CC3D5##7ac5cd"# "#535F7C"  #0063B2FF" #"#619CFF" "#9CC3D5FF"
      fam_lines_pos =  "grey80"##FCB42C" 
      fam_lines_neg = "grey30"##535F7C"
      ggplot(o, aes(x = year_, y = DPR)) + 
            geom_smooth(aes(color = slope_dir, fill = slope_dir), method = 'lm',  size = 0.5) +
            scale_color_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
            scale_fill_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  

            geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +

            new_scale_color() +   
            geom_point(aes(col = direct), alpha = 0.7, size = 1.1) + 
            scale_color_manual(name = "Directly\ncalculated\nestimates:",
                               values = c(direct, trans),
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

      ggsave(file="Outputs/Yearly_pop_points-lines-col2.png", dpi = 600, width = 14, height = 17, unit = "cm")
    # lines colored, points circles filled/empty
      direct = 19
      trans = 1
      fam_lines_pos =  "#FCB42C" 
      fam_lines_neg = "#535F7C"
      p_col = "#bc507d" #7dbc50"#grey50"##5eab2b"
      
      ggplot(o, aes(x = year_, y = DPR)) + 
            geom_smooth(aes(color = slope_dir, fill = slope_dir), method = 'lm',  size = 0.5) +
            scale_color_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
            scale_fill_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  

            geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +

            new_scale_color() +   
            geom_point(aes(shape = direct), size = 1.1, col = p_col) + 
            scale_shape_manual(name = "Directly\ncalculated\nestimates:",
                               values = c(direct, trans),
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

      ggsave(file="Outputs/Yearly_pop-points-circ_line-col.png", dpi = 600, width = 14, height = 17, unit = "cm")
    # lines colored, points shape
      direct = 19
      trans = 17
      fam_lines_pos =  "#FCB42C"  
      fam_lines_neg = "#535F7C" 
      p_col = "#bc507d" #7dbc50"#grey50"##5eab2b"
      
      ggplot(o, aes(x = year_, y = DPR)) + 
            geom_smooth(aes(color = slope_dir, fill = slope_dir), method = 'lm',  size = 0.5) +
            scale_color_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  
            scale_fill_manual(values=c(fam_lines_pos, fam_lines_neg), name = "Slope")+  

            geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +

            new_scale_color() +   
            geom_point(aes(shape = direct), size = 1.1, col = p_col) + 
            scale_shape_manual(name = "Directly\ncalculated\nestimates:",
                               values = c(direct, trans),
                               breaks = c("Yes", "No"))+

            facet_wrap(~labeler, ncol = 3, scale = "free_y", as.table = TRUE)+ 
            #
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

      ggsave(file="Outputs/Yearly_pop_points-shape_line-col.png", dpi = 600, width = 14, height = 17, unit = "cm")
  # with points colored
    direct = "#FCB42C" #"#5eab2b"#c3dbe5" #9CC3D5""#FCB42C"   #ffa500
    trans = "#535F7C"#"#7c9caa" #c3dbe5" #9CC3D5##7ac5cd"# "#535F7C"  #0063B2FF" #"#619CFF" "#9CC3D5FF"

    ggplot(o, aes(x = year_, y = DPR)) + 
          geom_smooth(method="lm", col = "black", lwd=0.5) +
          geom_vline(xintercept = 2000, lty=3, col = "grey80", lwd = 0.25) +
          geom_point(aes(col = direct), pch =16) + 
          facet_wrap(~labeler, ncol = 3, scale = "free_y", as.table = TRUE)+ 
          ylab("Daily predation rate") + xlab("Year") +
          scale_color_manual(name = "Directly\ncalculated\nestimates:",
                             values = c(direct, trans),
                             breaks = c("Yes", "No"))+
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
    ggsave(file="Outputs/Yearly_pop_points-col.png", dpi = 600, width = 14, height = 17, unit = "cm")
  # general
    ggplot(bbd, aes(x = year_, y = log(DPR+0.01))) + geom_point() + geom_smooth()
    ggplot(bbd[bbd$Belt!="South temperate",], aes(x = year_, y = log(DPR+0.01))) + geom_point() + geom_smooth() 
    ggplot(bbd[bbd$Belt!="South temperate",], aes(x = year_, y = log(DPR+0.01), col = Belt, fill = Belt)) + geom_point() + geom_smooth()
     

