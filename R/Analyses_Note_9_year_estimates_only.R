# TOOLS
  PNG = FALSE # print figures in PNG or not
  rm( list = ls() )	
  sapply(c('AICcmodavg','ape','arm','coxme','data.table','effects', 'ggnewscale','ggplot2','grid', 'lattice','magrittr','mgcv','multcomp','performance','phytools','plyr','RColorBrewer','readxl','writexl'), #XLConnect
      function(x) suppressPackageStartupMessages(require(x , character.only = TRUE, quietly = TRUE) ))
  source('R/Constants_Functions.R')

# ALL DATA
  # DATA preparation
    d = data.table(read.csv("Data/Rebuttal_data_with_DPR_by_yr.csv", h = T, sep=",",stringsAsFactors = FALSE))
    #dd = d[year!='all' & duplicated(paste(source_id,species, locality, year)), .(source_id, species, locality,year)] 
    #print(d[paste(source_id,species, locality, year)%in%dd[,paste(source_id,species, locality, year)],.(who, source_id,species, locality, Latitude, year, N.nests)], nrows =200)

    d[,pk := 1:nrow(d)]
    d[,Belt := as.factor(Belt)]
    d[,abs_lat := abs(Latitude)]
    d[,hemisphere :=as.factor(ifelse(Latitude > 0, "Northern", "Southern"))]
    d[,site := paste(Latitude, Longitude)]
    d = d[- which(is.na(year)|is.na(DPR))]
    #d[is.na(DPR)]
    d[,year_:=as.numeric(ifelse(year=='all', ((start_year+end_year)/2), d$year))]
    d[, N_nests_loc_year := sum(N.nests), by =list(species, locality, year)] # localities are larger or there are estimates for various habitats, so there are sometimes more estimates per locality in a given year
    
    #d[year!='all' & N_nests_loc_year!=N.nests,.(source_id,species, locality, Latitude, year, N.nests, N_nests_loc_year) ]
    d = d[N.nests>11] # uses more than 11 nests for fine scaled localities and the output is same as for d = d[d$N_nests_loc_year>11,] (except for Vanellus_vanellus source_id 118 where additional estimates based on 6,8, and 9 would be included) 
    #print(d[paste(source_id,species, locality, year)%in%dd[,paste(source_id,species, locality, year)], .(who, source_id,species, locality, year, N.nests,N_nests_loc_year)], nrows =200) 

    summary(factor(d$DPR_assumption))
    d[is.na(DPR_assumption), DPR_assumption:='our_est']

    b = data.table(d[which(d$year!='all'),])
    summary(factor(b$DPR_assumption))
    b[,fail_is_depr := ifelse(DPR_assumption == 'our_est', 'No','Yes')]

    b$year = as.numeric(b$year)
    b$n=1
    b$pop = paste(b$source_id,b$species,b$locality)
    b$direct = ifelse(b$DPRtrans == "YES", "No", "Yes")
    x = ddply(b,.(pop), summarise, n =length(unique(year)))
    summary(factor(x$n))
    nrow(x[x$n>5,])
    #x2 = ddply(b,.(pop), summarise, n =sum(n))
    #nrow(x2[x2$n>5,])
    #d[!pk%in%b$pk[b$species=='Vanellus_vanellus'] & species=='Vanellus_vanellus']

    o = data.table(b[b$pop%in%x$pop[x$n>5],])
    o[,  slope := lm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_dir := ifelse(slope >0, '+', '-')]
     o[,  slope_rlm_z := rlm(scale(DPR) ~ scale(year_))  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_rlm := rlm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_dir_rlm := ifelse(slope_rlm >0, '+', '-')]
    o[, labeler := factor(paste(round(Latitude), "|", species))]
    o[, labeler := factor(labeler, levels=rev(levels(labeler)))]

    oo =o[!duplicated(labeler)]
    summary(factor(oo$slope_dir))
    xtabs(~oo$slope_dir+oo$fail_is_depr+oo$direct) 
    summary(oo$slope_rlm_z)
    t.test(oo$slope_rlm[oo$Belt == 'Arctic'], oo$slope_rlm[oo$Belt != 'Arctic'])

  # TABLE
    xx = ddply(b,.(pop, locality), summarise, n =length(unique(year)))
    xy = ddply(xx,.(locality), summarise, n =sum(n))
    summary(factor(xy$n))

    bb= b[b$pop%in%x$pop[x$n>1],]
    bb$years_nr = x$n[match(bb$pop,x$pop)]
    xx = ddply(bb,.(pop, locality), summarise, n =length(unique(year)))
    xy = ddply(xx,.(locality), summarise, n =sum(n))
    summary(factor(xy$n))
   
    nrow(bb)
    length(unique(bb$pop))

    bd = b[b$DPRtrans=='NO',]
    bbd = bb[bb$DPRtrans=='NO',]

    nrow(b); length(unique(b$year_)); length(unique(b$site)); length(unique(b$species))
    nrow(bd); length(unique(bd$year_)); length(unique(bd$locality)); length(unique(bd$species))
    nrow(bb); length(unique(bb$year_)); length(unique(bb$locality)); length(unique(bb$species))
    nrow(bbd); length(unique(bbd$year_)); length(unique(bbd$locality)); length(unique(bbd$species))

    m1 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),b)
    m1b = lmer(log(DPR+0.01)~fail_is_depr+ log(N.nests)+year_ +(1|species) + (1|locality),b)
    m2 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bd)
    m2b = lmer(log(DPR+0.01)~fail_is_depr+log(N.nests)+year_ +(1|species) + (1|locality),bd)
    m3 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bb)
    m3b = lmer(log(DPR+0.01)~fail_is_depr+log(N.nests)+year_ +(1|species) + (1|locality),bb)
    m4 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bbd)
    m4b = lmer(log(DPR+0.01)~fail_is_depr+log(N.nests)+year_ +(1|species) + (1|locality),bbd)
        #summary(glht(m))
        #plot(allEffects(m)) 
    
    o1 = m_out(name = "all", model = m1, round_ = 5, nsim = 5000, aic = FALSE)  
    o2 = m_out(name = "good", model = m2, round_ = 3, nsim = 5000, aic = FALSE)  
    o3 = m_out(name = "all >1yr", model = m3, round_ = 3, nsim = 5000, aic = FALSE)  
    o4 = m_out(name = "good >1yr", model = m4, round_ = 3, nsim = 5000, aic = FALSE)  

    sname = 'Table_Note9_update_2020-12-09'
    write_xlsx(rbind(o1,o2,o3,o4), paste0("Outputs/",sname,'.xlsx'))

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

    ggsave(file="Outputs/rlm_Yearly_pop_points-lines-col-shapeEst.png", dpi = 600, width = 14, height = 17, unit = "cm")

# without years 2004,2005,2006,2007 from Lower Katanga alternative site.
  # DATA preparation
    d = data.table(read.csv("Data/Rebuttal_data_with_DPR_by_yr.csv", h = T, sep=",",stringsAsFactors = FALSE))
    d = d[!(source_id==208 & year%in%c(2004,2005,2006,2007))] # exludes data from Lower Katanga that were collected at a different field site.

    d[,pk := 1:nrow(d)]
    d[,Belt := as.factor(Belt)]
    d[,abs_lat := abs(Latitude)]
    d[,hemisphere :=as.factor(ifelse(Latitude > 0, "Northern", "Southern"))]
    d[,site := paste(Latitude, Longitude)]
    d = d[- which(is.na(year)|is.na(DPR))]
    #d[is.na(DPR)]
    d[,year_:=as.numeric(ifelse(year=='all', ((start_year+end_year)/2), d$year))]
    d[, N_nests_loc_year := sum(N.nests), by =list(species, locality, year)] # localities are larger or there are estimates for various habitats, so there are sometimes more estimates per locality in a given year
    
    #d[year!='all' & N_nests_loc_year!=N.nests,.(source_id,species, locality, Latitude, year, N.nests, N_nests_loc_year) ]
    d = d[N.nests>11] # uses more than 11 nests for fine scaled localities and the output is same as for d = d[d$N_nests_loc_year>11,] (except for Vanellus_vanellus source_id 118 where additional estimates based on 6,8, and 9 would be included) 
    #print(d[paste(source_id,species, locality, year)%in%dd[,paste(source_id,species, locality, year)], .(who, source_id,species, locality, year, N.nests,N_nests_loc_year)], nrows =200) 

    summary(factor(d$DPR_assumption))
    d[is.na(DPR_assumption), DPR_assumption:='our_est']

    b = data.table(d[which(d$year!='all'),])
    summary(factor(b$DPR_assumption))
    b[,fail_is_depr := ifelse(DPR_assumption == 'our_est', 'No','Yes')]

    b$year = as.numeric(b$year)
    b$n=1
    b$pop = paste(b$source_id,b$species,b$locality)
    b$direct = ifelse(b$DPRtrans == "YES", "No", "Yes")
    x = ddply(b,.(pop), summarise, n =length(unique(year)))
    summary(factor(x$n))
    nrow(x[x$n>5,])
    #x2 = ddply(b,.(pop), summarise, n =sum(n))
    #nrow(x2[x2$n>5,])
    #d[!pk%in%b$pk[b$species=='Vanellus_vanellus'] & species=='Vanellus_vanellus']

    o = data.table(b[b$pop%in%x$pop[x$n>5],])
    o[,  slope := lm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_dir := ifelse(slope >0, '+', '-')]
    o[,  slope_rlm := rlm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_dir_rlm := ifelse(slope_rlm >0, '+', '-')]
    o[, labeler := factor(paste(round(Latitude), "|", species))]
    o[, labeler := factor(labeler, levels=rev(levels(labeler)))]

    oo =o[!duplicated(labeler)]
    summary(factor(oo$slope_dir))
    xtabs(~oo$slope_dir+oo$fail_is_depr+oo$direct) 
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

    sname = 'Table_Note9_update_2020-12-09_withoutKatanga'
    write_xlsx(rbind(o1,o2,o3,o4), paste0("Outputs/",sname,'.xlsx'))
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

    ggsave(file="Outputs/rlm_Yearly_pop_points-lines-col_withoutAltKatanga.png", dpi = 600, width = 14, height = 17, unit = "cm")

# only DPR based on predation (not including all failures)
  # DATA preparation
    d = data.table(read.csv("Data/Rebuttal_data_with_DPR_by_yr.csv", h = T, sep=",",stringsAsFactors = FALSE))

    d[,pk := 1:nrow(d)]
    d[,Belt := as.factor(Belt)]
    d[,abs_lat := abs(Latitude)]
    d[,hemisphere :=as.factor(ifelse(Latitude > 0, "Northern", "Southern"))]
    d[,site := paste(Latitude, Longitude)]
    d = d[- which(is.na(year)|is.na(DPR))]
    #d[is.na(DPR)]
    d[,year_:=as.numeric(ifelse(year=='all', ((start_year+end_year)/2), d$year))]
    d[, N_nests_loc_year := sum(N.nests), by =list(species, locality, year)] # localities are larger or there are estimates for various habitats, so there are sometimes more estimates per locality in a given year
    
    #d[year!='all' & N_nests_loc_year!=N.nests,.(source_id,species, locality, Latitude, year, N.nests, N_nests_loc_year) ]
    d = d[N.nests>11] # uses more than 11 nests for fine scaled localities and the output is same as for d = d[d$N_nests_loc_year>11,] (except for Vanellus_vanellus source_id 118 where additional estimates based on 6,8, and 9 would be included) 
    #print(d[paste(source_id,species, locality, year)%in%dd[,paste(source_id,species, locality, year)], .(who, source_id,species, locality, year, N.nests,N_nests_loc_year)], nrows =200) 

    summary(factor(d$DPR_assumption))
    d = d[is.na(DPR_assumption)]

    b = data.table(d[which(d$year!='all'),])
    
    b$year = as.numeric(b$year)
    b$n=1
    b$pop = paste(b$source_id,b$species,b$locality)
    b$direct = ifelse(b$DPRtrans == "YES", "No", "Yes")
    x = ddply(b,.(pop), summarise, n =length(unique(year)))
    summary(factor(x$n))
    nrow(x[x$n>5,])
    #x2 = ddply(b,.(pop), summarise, n =sum(n))
    #nrow(x2[x2$n>5,])
    #d[!pk%in%b$pk[b$species=='Vanellus_vanellus'] & species=='Vanellus_vanellus']

    o = data.table(b[b$pop%in%x$pop[x$n>5],])
    o[,  slope := lm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_dir := ifelse(slope >0, '+', '-')]
    o[,  slope_rlm_z := rlm(scale(DPR) ~ scale(year_))  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_rlm := rlm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_dir_rlm := ifelse(slope_rlm >0, '+', '-')]
    o[, labeler := factor(paste(round(Latitude), "|", species))]
    o[, labeler := factor(labeler, levels=rev(levels(labeler)))]

   oo =o[!duplicated(labeler)]
    summary(factor(oo$slope_dir))
    xtabs(~oo$slope_dir+oo$direct) 
    summary(oo$slope_rlm_z)
    t.test(oo$slope_rlm[oo$Belt == 'Arctic'], oo$slope_rlm[oo$Belt != 'Arctic'])

  # TABLE
    bb= b[b$pop%in%x$pop[x$n>1],]
    bb$years_nr = x$n[match(bb$pop,x$pop)]
    xx = ddply(bb,.(pop, locality), summarise, n =length(unique(year)))
    xy = ddply(xx,.(locality), summarise, n =sum(n))
    summary(factor(xy$n))
   
    nrow(bb)
    length(unique(bb$pop))

    bd = b[b$DPRtrans=='NO',]
    bbd = bb[bb$DPRtrans=='NO',]
    xx = ddply(bbd,.(pop, locality), summarise, n =length(unique(year)))
    xy = ddply(xx,.(locality), summarise, n =sum(n))
    summary(factor(xy$n))
 
    nrow(b); length(unique(b$year_)); length(unique(b$site)); length(unique(b$species))
    nrow(bd); length(unique(bd$year_)); length(unique(bd$locality)); length(unique(bd$species))
    nrow(bb); length(unique(bb$year_)); length(unique(bb$locality)); length(unique(bb$species))
    nrow(bbd); length(unique(bbd$year_)); length(unique(bbd$locality)); length(unique(bbd$species))

    m1 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),b)
      m1 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),b)
    summary(glht(m1))
    m2 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bd)
      m2 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),bd)
     summary(glht(m2))
    m3 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bb)
      m3 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),bb)
     summary(glht(m3))
    m4 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bbd)
      m4 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),bbd)
        summary(glht(m4))
        #plot(allEffects(m)) 
    
    o1 = m_out(name = "all", model = m1, round_ = 3, nsim = 5000, aic = FALSE)  
    o2 = m_out(name = "good", model = m2, round_ = 3, nsim = 5000, aic = FALSE)  
    o3 = m_out(name = "all >1yr", model = m3, round_ = 3, nsim = 5000, aic = FALSE)  
    o4 = m_out(name = "good >1yr", model = m4, round_ = 3, nsim = 5000, aic = FALSE)  

    sname = 'Table_Note9_update_2020-12-09_onlyPred_noOtherFailures'
    write_xlsx(rbind(o1,o2,o3,o4), paste0("Outputs/",sname,'.xlsx'))

    #ggplot(bbd, aes(x = year_, y = DPR))+geom_smooth() +geom_point(alpha = 0.3)
    #ggplot(bbd, aes(x = year_, y = log(DPR+0.01)))+geom_smooth() +geom_point(alpha = 0.3)
    #bbd[Belt == 'Arctic', arctic :='Yes']
    #bbd[Belt != 'Arctic', arctic :='No']
    #m4 = lmer(log(DPR+0.01)~log(N.nests)+arctic*scale(year_) +(1|species) + (scale(year_)|locality),bbd)
    #m4 = lmer(log(DPR+0.01)~log(N.nests)+arctic*poly(year_,2) +(1|species) + (poly(year_,2)|locality),bbd)
  # TABLE with slopes
    bb= b[b$pop%in%x$pop[x$n>1],]
    bb$years_nr = x$n[match(bb$pop,x$pop)]
    xx = ddply(bb,.(pop, locality), summarise, n =length(unique(year)))
    xy = ddply(xx,.(locality), summarise, n =sum(n))
    summary(factor(xy$n))
   
    nrow(bb)
    length(unique(bb$pop))

    bd = b[b$DPRtrans=='NO',]
    bbd = bb[bb$DPRtrans=='NO',]
    xx = ddply(bbd,.(pop, locality), summarise, n =length(unique(year)))
    xy = ddply(xx,.(locality), summarise, n =sum(n))
    summary(factor(xy$n))
 
    nrow(b); length(unique(b$year_)); length(unique(b$site)); length(unique(b$species))
    nrow(bd); length(unique(bd$year_)); length(unique(bd$locality)); length(unique(bd$species))
    nrow(bb); length(unique(bb$year_)); length(unique(bb$locality)); length(unique(bb$species))
    nrow(bbd); length(unique(bbd$year_)); length(unique(bbd$locality)); length(unique(bbd$species))

    #m1 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),b)
    m1 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),b)
    summary(glht(m1))
    #m2 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bd)
    m2 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),bd)
     summary(glht(m2))
    #m3 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bb)
    m3 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),bb)
     summary(glht(m3))
    #m4 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bbd)
    m4 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),bbd)
        summary(glht(m4))
        #plot(allEffects(m)) 
    
    o1 = m_out(name = "all", model = m1, round_ = 3, nsim = 5000, aic = FALSE)  
    o2 = m_out(name = "good", model = m2, round_ = 3, nsim = 5000, aic = FALSE)  
    o3 = m_out(name = "all >1yr", model = m3, round_ = 3, nsim = 5000, aic = FALSE)  
    o4 = m_out(name = "good >1yr", model = m4, round_ = 3, nsim = 5000, aic = FALSE)  

    sname = 'Table_Note9_update_2020-12-09_onlyPred_noOtherFailures_slopes'
    write_xlsx(rbind(o1,o2,o3,o4), paste0("Outputs/",sname,'.xlsx'))

    #ggplot(bbd, aes(x = year_, y = DPR))+geom_smooth() +geom_point(alpha = 0.3)
    #ggplot(bbd, aes(x = year_, y = log(DPR+0.01)))+geom_smooth() +geom_point(alpha = 0.3)
    #bbd[Belt == 'Arctic', arctic :='Yes']
    #bbd[Belt != 'Arctic', arctic :='No']
    #m4 = lmer(log(DPR+0.01)~log(N.nests)+arctic*scale(year_) +(1|species) + (scale(year_)|locality),bbd)
    #m4 = lmer(log(DPR+0.01)~log(N.nests)+arctic*poly(year_,2) +(1|species) + (poly(year_,2)|locality),bbd)
  
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
          geom_point(aes(col = direct), size = 1.1, shape = 16) + 
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

    ggsave(file="Outputs/rlm_Yearly_pop_points-lines-col-shapeEst_BestData.png", dpi = 600, width = 14, height = 17, unit = "cm")
# only DPR based on predation (not including all failures) and only DIRECTLY calculated data
  # DATA preparation
    d = data.table(read.csv("Data/Rebuttal_data_with_DPR_by_yr.csv", h = T, sep=",",stringsAsFactors = FALSE))

    d[,pk := 1:nrow(d)]
    d[,Belt := as.factor(Belt)]
    d[,abs_lat := abs(Latitude)]
    d[,hemisphere :=as.factor(ifelse(Latitude > 0, "Northern", "Southern"))]
    d[,site := paste(Latitude, Longitude)]
    d = d[- which(is.na(year)|is.na(DPR))]
    #d[is.na(DPR)]
    d[,year_:=as.numeric(ifelse(year=='all', ((start_year+end_year)/2), d$year))]
    d[, N_nests_loc_year := sum(N.nests), by =list(species, locality, year)] # localities are larger or there are estimates for various habitats, so there are sometimes more estimates per locality in a given year
    
    #d[year!='all' & N_nests_loc_year!=N.nests,.(source_id,species, locality, Latitude, year, N.nests, N_nests_loc_year) ]
    d = d[N.nests>11] # uses more than 11 nests for fine scaled localities and the output is same as for d = d[d$N_nests_loc_year>11,] (except for Vanellus_vanellus source_id 118 where additional estimates based on 6,8, and 9 would be included) 
    #print(d[paste(source_id,species, locality, year)%in%dd[,paste(source_id,species, locality, year)], .(who, source_id,species, locality, year, N.nests,N_nests_loc_year)], nrows =200) 

    summary(factor(d$DPR_assumption))
    d = d[is.na(DPR_assumption)]

    b = data.table(d[which(d$year!='all'),])
    
    b$year = as.numeric(b$year)
    b$n=1
    b$pop = paste(b$source_id,b$species,b$locality)
    b$direct = ifelse(b$DPRtrans == "YES", "No", "Yes")
    b = b[direct == 'Yes']
    x = ddply(b,.(pop), summarise, n =length(unique(year)))
    summary(factor(x$n))
    nrow(x[x$n>5,])
    #x2 = ddply(b,.(pop), summarise, n =sum(n))
    #nrow(x2[x2$n>5,])
    #d[!pk%in%b$pk[b$species=='Vanellus_vanellus'] & species=='Vanellus_vanellus']

    o = data.table(b[b$pop%in%x$pop[x$n>5],])
    o[,  slope := lm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_dir := ifelse(slope >0, '+', '-')]
    o[,  slope_rlm_z := rlm(scale(DPR) ~ scale(year_))  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_rlm := rlm(DPR ~ year_)  %>% coef  %>% extract(2), by = pop] 
    o[,  slope_dir_rlm := ifelse(slope_rlm >0, '+', '-')]
    o[, labeler := factor(paste(round(Latitude), "|", species))]
    o[, labeler := factor(labeler, levels=rev(levels(labeler)))]

   oo =o[!duplicated(labeler)]
    summary(factor(oo$slope_dir))
    xtabs(~oo$slope_dir+oo$direct) 
    summary(oo$slope_rlm_z)
    summary(factor(oo$Belt)) # only two non-arctic sites
    t.test(oo$slope_rlm[oo$Belt == 'Arctic'], oo$slope_rlm[oo$Belt != 'Arctic'])

  # TABLE
    bb= b[b$pop%in%x$pop[x$n>1],]
    bb$years_nr = x$n[match(bb$pop,x$pop)]
    xx = ddply(bb,.(pop, locality), summarise, n =length(unique(year)))
    xy = ddply(xx,.(locality), summarise, n =sum(n))
    summary(factor(xy$n))
   
    nrow(bb)
    length(unique(bb$pop))

    bd = b[b$DPRtrans=='NO',]
    bbd = bb[bb$DPRtrans=='NO',]
    xx = ddply(bbd,.(pop, locality), summarise, n =length(unique(year)))
    xy = ddply(xx,.(locality), summarise, n =sum(n))
    summary(factor(xy$n))
 
    nrow(b); length(unique(b$year_)); length(unique(b$site)); length(unique(b$species))
    nrow(bd); length(unique(bd$year_)); length(unique(bd$locality)); length(unique(bd$species))
    nrow(bb); length(unique(bb$year_)); length(unique(bb$locality)); length(unique(bb$species))
    nrow(bbd); length(unique(bbd$year_)); length(unique(bbd$locality)); length(unique(bbd$species))

    m1 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),b)
      m1 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),b)
    summary(glht(m1))
    m2 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bd)
      m2 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),bd)
     summary(glht(m2))
    m3 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bb)
      m3 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),bb)
     summary(glht(m3))
    m4 = lmer(log(DPR+0.01)~log(N.nests)+year_ +(1|species) + (1|locality),bbd)
      m4 = lmer(log(DPR+0.01)~log(N.nests)+scale(year_) +(1|species)+ (scale(year_)|locality),bbd)
        summary(glht(m4))
        #plot(allEffects(m)) 
    
    o1 = m_out(name = "all", model = m1, round_ = 3, nsim = 5000, aic = FALSE)  
    o2 = m_out(name = "good", model = m2, round_ = 3, nsim = 5000, aic = FALSE)  
    o3 = m_out(name = "all >1yr", model = m3, round_ = 3, nsim = 5000, aic = FALSE)  
    o4 = m_out(name = "good >1yr", model = m4, round_ = 3, nsim = 5000, aic = FALSE)  

    sname = 'Table_Note9_update_2020-12-09_onlyPred_noOtherFailures'
    write_xlsx(rbind(o1,o2,o3,o4), paste0("Outputs/",sname,'.xlsx'))

    #ggplot(bbd, aes(x = year_, y = DPR))+geom_smooth() +geom_point(alpha = 0.3)
    #ggplot(bbd, aes(x = year_, y = log(DPR+0.01)))+geom_smooth() +geom_point(alpha = 0.3)
    #bbd[Belt == 'Arctic', arctic :='Yes']
    #bbd[Belt != 'Arctic', arctic :='No']
    #m4 = lmer(log(DPR+0.01)~log(N.nests)+arctic*scale(year_) +(1|species) + (scale(year_)|locality),bbd)
    #m4 = lmer(log(DPR+0.01)~log(N.nests)+arctic*poly(year_,2) +(1|species) + (poly(year_,2)|locality),bbd)
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
          geom_point(aes(col = direct), size = 1.1, shape = 16) + 
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

    ggsave(file="Outputs/rlm_Yearly_pop_points-lines-col-shapeEst_BestData&DirectOnly.png", dpi = 600, width = 14, height = 17, unit = "cm")


# PLOT per site
  # DATA preparation
    d = data.table(read.csv("Data/Rebuttal_data_with_DPR_by_yr.csv", h = T, sep=",",stringsAsFactors = FALSE))

    d[,pk := 1:nrow(d)]
    d[,Belt := as.factor(Belt)]
    d[,abs_lat := abs(Latitude)]
    d[,hemisphere :=as.factor(ifelse(Latitude > 0, "Northern", "Southern"))]
    d = d[- which(is.na(year)|is.na(DPR))]
    #d[is.na(DPR)]
    d[,year_:=as.numeric(ifelse(year=='all', ((start_year+end_year)/2), d$year))]
    d = d[N.nests>11] # uses more than 11 nests for fine scaled localities and the output is same as for d = d[d$N_nests_loc_year>11,] (except for Vanellus_vanellus source_id 118 where additional estimates based on 6,8, and 9 would be included) 

    summary(factor(d$DPR_assumption))
    d[is.na(DPR_assumption), DPR_assumption:='our_est']

    b = data.table(d[which(d$year!='all'),])
    summary(factor(b$DPR_assumption))
    b[,fail_is_depr := ifelse(DPR_assumption == 'our_est', 'No','Yes')]

    b$year = as.numeric(b$year)
    b$n=1
    b[,labeler := paste(round(Latitude,2), "|", round(Longitude,2))]
    b$direct = ifelse(b$DPRtrans == "YES", "No", "Yes")
    x = ddply(b,.(labeler), summarise, n =length(unique(year)))
    summary(factor(x$n))
    nrow(x[x$n>5,])
    
    o = data.table(b[b$labeler%in%x$labeler[x$n>5],])
    o[,  slope := lm(DPR ~ year_)  %>% coef  %>% extract(2), by = labeler] 
    o[,  slope_dir := ifelse(slope >0, '+', '-')]
    o[,  slope_rlm := rlm(DPR ~ year_)  %>% coef  %>% extract(2), by = labeler] 
    o[,  slope_dir_rlm := ifelse(slope_rlm >0, '+', '-')]
    o[, labeler := factor(labeler)]
    o[, labeler := factor(labeler, levels=rev(levels(labeler)))]
  # PLOT all
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

    ggsave(file="Outputs/rlm_Yearly_site.png", dpi = 600, width = 14, height = 17, unit = "cm")
  # PLOT means per year
    o = data.table(b[b$labeler%in%x$labeler[x$n>5],]) 
    o = data.table(ddply(o,.(labeler, year_,direct,fail_is_depr ),summarise,  DPR = mean(DPR)))
    o[,  slope := lm(DPR ~ year_)  %>% coef  %>% extract(2), by = labeler] 
    o[,  slope_dir := ifelse(slope >0, '+', '-')]
    o[,  slope_rlm := rlm(DPR ~ year_)  %>% coef  %>% extract(2), by = labeler] 
    o[,  slope_dir_rlm := ifelse(slope_rlm >0, '+', '-')]
    o[, labeler := factor(labeler)]
    o[, labeler := factor(labeler, levels=rev(levels(labeler)))]

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
          ylab("Mean daily predation rate") + xlab("Year") +
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

    ggsave(file="Outputs/rlm_Yearly_site_means.png", dpi = 600, width = 14, height = 17, unit = "cm")

##################

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
     

