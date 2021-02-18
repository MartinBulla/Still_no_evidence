# ===========================================================================
# Supporting information for "Bulla, Valcu & Kempenaers. (2021) "Still no 
# evidence for disruption of global patterns of nest predation in shorebirds" 
# Contributor: Martin Bulla
# 📍 script runs relative to the project's root directory, uses Kubelka et 
# al's data and its manipulation, and generates Figure R2 within part
# (1) Multiple sampling and pseudo-replication
# at places Kubelka et al have used log(dpr+0.01), but the dpr variable 
# contains “+0.01” already. Thus, they have done log(dpr+0.02) and for 
# consistency we do the same.
# ===========================================================================

# PREPARATION
  rm(list = ls())

  # packages
  library(arm) 
  library(MuMIn) 
  library(lme4)
  library(lmerTest)
  library(AICcmodavg)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(multcomp)

  # Calculate standard error
  stderr <- function(x) sd(x) / sqrt(length(x))

  # data
    datapred <- read.csv("Data/DATApopulations.csv", sep =";", h = T)
    (datapred$hemis <-  datapred$Latitude > 0 )
    (datapred$hemisphere <- ifelse(datapred$Latitude > 0, 'North','South' ))
    (datapred$isarctic <- datapred$Latitude > 60)
    (datapred$decade <- trunc(datapred$mean_year/10))
    (datapred$pop_ID <-paste(datapred$Longitude, datapred$Latitude))
    datapred$absLat = abs(datapred$Latitude)

    summary(datapred$mean_year[datapred$Belt == 'Arctic']) 
    summary(datapred$mean_year[datapred$Belt == 'North temperate']) 
    
# Figure R2 top - Kubelka et al's C, D Arctic versus North Temperate

  # prepare data summary 
    # means - as in Kubelka et al
      (dataRed <- datapred[c(which(datapred$Belt == "North temperate"), which(datapred$Belt == "Arctic")),])
       (datasummary <- dataRed %>% group_by( isarctic, decade) %>% 
        summarize( meandpr = exp( mean(log(DPR)  )),
                   se = stderr( DPR ),
                   plusse =  exp( mean(log(DPR)  ) + stderr( log(DPR) )   ),
                   minusse =  exp( mean(log(DPR)  ) -  stderr( log(DPR) ) ),
                   sel  = stderr( log(DPR) ),
                    n = length(DPR) ) )
    
    # prepare weighted means
      (datasummaryW <- dataRed %>% group_by( isarctic, decade) %>% 
        summarize( meandpr = exp( weighted.mean(log(DPR), N_nests )),
                   se = stderr( DPR ),
                   plusse =  exp( exp( weighted.mean(log(DPR), N_nests )) + stderr( log(DPR) )   ),
                   minusse =  exp(exp( weighted.mean(log(DPR), N_nests )) -  stderr( log(DPR) ) ),
                   sel  = stderr( log(DPR) ),
                   n = length(DPR) ) )  
    # prepare weighted dataset per site and then average 
      (datasum <- dataRed %>% group_by( isarctic, decade, pop_ID) %>% 
        summarize( meandpr = exp( weighted.mean(log(DPR), N_nests )),
                   n = length(DPR) ) )

      (datasummaryWW <- datasum %>% group_by( isarctic, decade) %>% 
        summarize( meandpr = exp( weighted.mean(log(meandpr), n )),
                   n = length(isarctic) ) )

  # explore differences in various summary methods  
    p1 =ggplot() + geom_point(aes(x = datasummaryW$meandpr, y = datasummary$meandpr)) + geom_abline(yintercpt = 0, slope = 1)
    p2 = ggplot() + geom_point(aes(x = datasummaryWW$meandpr, y = datasummary$meandpr)) + geom_abline(yintercpt = 0, slope = 1)
    p3 =ggplot() + geom_point(aes(x = datasummaryWW$meandpr, y = datasummaryW$meandpr)) + geom_abline(yintercpt = 0, slope = 1)
    grid.arrange(p1,p2,p3)

    f1 = ggplot( datasummary, aes( x = decade*10+5, y = meandpr ) ) + aes( colour = isarctic) +
      labs(x = "Year", y = "Mean DPR") +
      ylim(0, 0.07)+ 
      geom_point( aes(size =n), position = position_dodge(width = 0.25)) +
      scale_color_manual(breaks = c( "FALSE", "TRUE"), values = c("darkolivegreen4", "deeppink" ))+
      theme_bw() + 
      theme( panel.border = element_blank(), 
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             #legend.position="none",
             axis.line = element_line(colour = 'black', size = 0.25),
             axis.ticks = element_line(colour = "black", size = 0.25),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 8),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size = 10) )

      f2 = ggplot( datasummaryW, aes( x = decade*10+5, y = meandpr ) ) + aes( colour = isarctic) +
      labs(x = "Year", y = "Mean DPR") +
      ylim(0, 0.07)+ 
      geom_point( aes(size =n), position = position_dodge(width = 0.25)) +
      scale_color_manual(breaks = c( "FALSE", "TRUE"), values = c("darkolivegreen4", "deeppink" ), guide = FALSE )+
      theme_bw() + 
      theme( panel.border = element_blank(), 
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             #legend.position="none",
             axis.line = element_line(colour = 'black', size = 0.25),
             axis.ticks = element_line(colour = "black", size = 0.25),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 8),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size = 10) )

     f3 = ggplot( datasummaryWW, aes( x = decade*10+5, y = meandpr ) ) + aes( colour = isarctic) +
      labs(x = "Year", y = "Mean DPR") +
      ylim(0, 0.07)+ 
      geom_point( aes(size =n), position = position_dodge(width = 0.25)) +
      scale_color_manual(breaks = c( "FALSE", "TRUE"), values = c("darkolivegreen4", "deeppink" ), guide = FALSE )+
      theme_bw() + 
      theme( panel.border = element_blank(), 
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             #legend.position="none",
             axis.line = element_line(colour = 'black', size = 0.25),
             axis.ticks = element_line(colour = "black", size = 0.25),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 8),
             axis.title.x=element_blank(),
             axis.title.y=element_text(size = 10) )

      grid.arrange(f1,f2,f3)
      ggsave('Outputs/respEmail_Fig1C_N_iweighted_means.png',f2, width = 3.5, height = 2)
  
  # plot Kubelka et al's original Figure C 
    model.lm1 <- lm( log(DPR + 0.01) ~ poly(mean_year, 1) + Latitude  , data = datapred) 
    years <- seq(1944, 2016, by = 1)

    (newarctic <- data.frame( mean_year = years, Latitude = rep(69, 73) ))
    (newtempN  <- data.frame( mean_year = years, Latitude = rep(50.78, 73) ))

    predarctic  <- exp( predict(model.lm1, newarctic) ) - 0.01
    predtempN <- exp(  predict(model.lm1, newtempN) ) - 0.01

    preda <- data.frame( years = years, predarctic = predarctic, isarctic = TRUE )
    predo <- data.frame( years = years, predtempN = predtempN, isarctic = FALSE )

    figure <- ggplot( datasummary, aes( x = decade*10+5, y = meandpr ) ) + aes( colour = isarctic) +
      labs(x = "Year", y = "Mean DPR") +
      ylim(0, 0.07)+ 
      geom_errorbar( aes(ymin = minusse, ymax = plusse ), width = 0.1, position = position_dodge(width = 0.25) )+
      geom_point( position = position_dodge(width = 0.25)) +
      theme_bw() + 
      theme( panel.border = element_blank(), 
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             #legend.position="none",
             axis.line = element_line(colour = 'black', size = 0.25),
             axis.ticks = element_line(colour = "black", size = 0.25),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 8),
             #axis.title.x=element_blank(),
             axis.title.y=element_text(size = 10) )
    figure <- figure + scale_color_manual(breaks = c( "FALSE", "TRUE"), values = c("darkolivegreen4", "deeppink" ), guide = FALSE)

    figure <- figure + geom_line( aes(x = years, y = predarctic), preda  ) + geom_line( aes(x = years, y = predtempN), predo  )

    figure = figure + geom_text(x=1953, y=0.0625, label="Arctic", size = 3, col ='grey30', adj = 0) + geom_point(x = 1950, y = 0.0625,col = "deeppink")
    figure = figure + geom_text(x=1953, y=0.0565, label="North temperate", col = "grey30", size = 3, adj = 0)+ geom_point(x = 1950, y = 0.0565,col = "darkolivegreen4")

    ggsave('Outputs/respEmail_Fig1C_original.png', figure, , width = 3.5, height = 2)

    f_title = figure + ggtitle(label ="Figure 1C as in Kubelka et al. 2019b", subtitle = "lines based on\nlm(log(DPR+0.01)~poly(mean_year,1)+Latitude") + 
        theme(plot.title = element_text(size=9),
          plot.subtitle = element_text(size=8))
    ggsave('Outputs/Figure_R2_topLeft_Fig1C_original.png', f_title, width = 3.2, height = 2.75)

  # plot predictions - controlled for population and using interaction with linear latitude (not absolute)
    # prepare predictions
      m <- lmer( log(DPR + 0.01) ~ scale(mean_year)*scale(Latitude) + (1|pop_ID) , data = datapred)  
      summary(glht(m))
      m <- lmer( log(DPR + 0.01) ~ mean_year*Latitude + (1|pop_ID) , data = datapred) 
      nsim <- 5000
      bsim <- sim(m, n.sim=nsim) 
      v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
              apply(bsim@fixef, 2, quantile, prob=c(0.025,0.5,0.975))
          # values to predict for   
          newD=data.frame(mean_year = seq(1968,2015, length.out=300), 
                          Latitude = 69)
           
          # exactly the model which was used has to be specified here 
            X <- model.matrix(~ mean_year*Latitude,data=newD)  
                      
          # calculate predicted values and creditability intervals
            newD$pred <- exp(X%*%v)-0.01 # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(i in 1:nsim) predmatrix[,i] <- exp(X%*%bsim@fixef[i,])-0.01
                newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
            arc=newD   

          # values to predict for   
          newD=data.frame(mean_year = seq(1944,2013, length.out=300), 
                          Latitude = 50.78)
           
          # exactly the model which was used has to be specified here 
            X <- model.matrix(~ mean_year*Latitude,data=newD)  
                      
          # calculate predicted values and creditability intervals
            newD$pred <- exp(X%*%v)-0.01 # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(i in 1:nsim) predmatrix[,i] <- exp(X%*%bsim@fixef[i,])-0.01
                newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
            NTemp=newD  

        yearsA <- arc$mean_year
        yearsT <- NTemp$mean_year

        predarctic  <- arc$pred
        predtempN <- NTemp$pred
        predarcticUPR = arc$upr
        predarcticLWR = arc$lwr
        predtempNUPR = NTemp$upr
        predtempNLWR  = NTemp$lwr

      preda <- data.frame( years = yearsA, predarctic = predarctic,isarctic = TRUE )
      predaLwr <- data.frame( years = yearsA, predarctic = predarcticLWR,isarctic = TRUE )
      predaUpr <- data.frame( years = yearsA, predarctic = predarcticUPR,isarctic = TRUE )
      predo <- data.frame( years = yearsT, predtempN = predtempN, isarctic = FALSE )
      predoLwr <- data.frame( years = yearsT, predtempN = predtempNUPR, isarctic = FALSE )
      predoUpr <- data.frame( years = yearsT, predtempN = predtempNLWR, isarctic = FALSE )
    # plot
      figurebLinLat <- ggplot( datasummary, aes( x = decade*10+5, y = meandpr ) ) + aes( colour = isarctic) +
        labs(x = "Year", y = "Mean DPR") +
        ylim(0, 0.07)+ 
        geom_errorbar( aes(ymin = minusse, ymax = plusse ), width = 0.1, position = position_dodge(width = 0.25) )+
        geom_point( aes(size =n), position = position_dodge(width = 0.25)) +
        guides(size=guide_legend(title=expression(N[estimates])))+
        theme_bw() + 
        theme( panel.border = element_blank(), 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               #legend.position="none",
               axis.line = element_line(colour = 'black', size = 0.25),
               axis.ticks = element_line(colour = "black", size = 0.25),
               axis.text.x = element_text(size = 10),
               axis.text.y = element_text(size = 8),
               #axis.title.x=element_blank(),
               axis.title.y=element_text(size = 10) )
      figurebLinLat <- figurebLinLat + scale_color_manual(breaks = c( "FALSE", "TRUE"), values = c("darkolivegreen4", "deeppink" ), guide = FALSE)
      figurebLinLat <- figurebLinLat + 
        geom_line( aes(x = years, y = predarctic), preda  ) + 
        geom_line( aes(x = years, y = predarctic), predaLwr, lty = 3) + 
        geom_line( aes(x = years, y = predarctic), predaUpr, lty = 3) + 

        geom_line( aes(x = years, y = predtempN), predo  ) +
        geom_line( aes(x = years, y = predtempN), predoLwr, lty = 3) + 
        geom_line( aes(x = years, y = predtempN), predoUpr, lty = 3)
      ggsave('Outputs/respEmail_Fig1C_N_interaction_siteControl_latNotabs.png',figurebLinLat, width = 3.2, height = 2.75) 

      fC_title_lat = figurebLinLat + ggtitle(label ="Test for interaction & control for site", subtitle = "lines based on\nlm(log(DPR+0.01)~mean_year*Latitude+(1|site)") + 
        theme(plot.title = element_text(size=9),
          plot.subtitle = element_text(size=8),
          legend.position="none")
      ggsave('Outputs/Figure_R2_topRight_Fig1C_N_interaction_siteControl_latNotabs.png', fC_title_lat, , width = 3.2, height = 2.75)  
   
  # not used - plot predictions - controlled for population and using interactions with absolute latitude 
    # prepare predictions
      m <- lmer( log(DPR + 0.01) ~ scale(mean_year)*scale(abs(Latitude)) + (1|pop_ID) , data = datapred)  
      summary(glht(m))
      m <- lmer( log(DPR + 0.01) ~ mean_year*absLat + (1|pop_ID) , data = datapred) 
      nsim <- 5000
      bsim <- sim(m, n.sim=nsim) 
      v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
              apply(bsim@fixef, 2, quantile, prob=c(0.025,0.5,0.975))
          # values to predict for   
          newD=data.frame(mean_year = seq(1968,2015, length.out=300), 
                          absLat = 69)
           
          # exactly the model which was used has to be specified here 
            X <- model.matrix(~ mean_year*absLat,data=newD)  
                      
          # calculate predicted values and creditability intervals
            newD$pred <- exp(X%*%v)-0.01 # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(i in 1:nsim) predmatrix[,i] <- exp(X%*%bsim@fixef[i,])-0.01
                newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
            arc=newD   

          # values to predict for   
          newD=data.frame(mean_year = seq(1944,2013, length.out=300), 
                          absLat = 50.78)
           
          # exactly the model which was used has to be specified here 
            X <- model.matrix(~ mean_year*absLat,data=newD)  
                      
          # calculate predicted values and creditability intervals
            newD$pred <- exp(X%*%v)-0.01 # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(i in 1:nsim) predmatrix[,i] <- exp(X%*%bsim@fixef[i,])-0.01
                newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
            NTemp=newD  

        yearsA <- arc$mean_year
        yearsT <- NTemp$mean_year

        predarctic  <- arc$pred
        predtempN <- NTemp$pred
        predarcticUPR = arc$upr
        predarcticLWR = arc$lwr
        predtempNUPR = NTemp$upr
        predtempNLWR  = NTemp$lwr

      preda <- data.frame( years = yearsA, predarctic = predarctic,isarctic = TRUE )
      predaLwr <- data.frame( years = yearsA, predarctic = predarcticLWR,isarctic = TRUE )
      predaUpr <- data.frame( years = yearsA, predarctic = predarcticUPR,isarctic = TRUE )
      predo <- data.frame( years = yearsT, predtempN = predtempN, isarctic = FALSE )
      predoLwr <- data.frame( years = yearsT, predtempN = predtempNUPR, isarctic = FALSE )
      predoUpr <- data.frame( years = yearsT, predtempN = predtempNLWR, isarctic = FALSE )
    # plot
      figureb <- ggplot( datasummary, aes( x = decade*10+5, y = meandpr ) ) + aes( colour = isarctic) +
        labs(x = "Year", y = "Mean DPR") +
        ylim(0, 0.07)+ 
        geom_errorbar( aes(ymin = minusse, ymax = plusse ), width = 0.1, position = position_dodge(width = 0.25) )+
        geom_point( aes(size =n), position = position_dodge(width = 0.25)) +
        guides(size=guide_legend(title=expression(N[estimates])))+
        theme_bw() + 
        theme( panel.border = element_blank(), 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               #legend.position="none",
               axis.line = element_line(colour = 'black', size = 0.25),
               axis.ticks = element_line(colour = "black", size = 0.25),
               axis.text.x = element_text(size = 10),
               axis.text.y = element_text(size = 8),
               #axis.title.x=element_blank(),
               axis.title.y=element_text(size = 10) )
      figureb <- figureb + scale_color_manual(breaks = c( "FALSE", "TRUE"), values = c("darkolivegreen4", "deeppink" ), guide = FALSE)
  
      figureb <- figureb + 
          geom_line( aes(x = years, y = predarctic), preda  ) + 
          geom_line( aes(x = years, y = predarctic), predaLwr, lty = 3) + 
          geom_line( aes(x = years, y = predarctic), predaUpr, lty = 3) + 

          geom_line( aes(x = years, y = predtempN), predo  ) +
          geom_line( aes(x = years, y = predtempN), predoLwr, lty = 3) + 
          geom_line( aes(x = years, y = predtempN), predoUpr, lty = 3)

      ggsave('Outputs/respEmail_Fig1C_N_interaction_siteControl_absLat.png',figureb, width = 3.2, height = 2)

      fC_title = figureb + ggtitle(label ="Test for interaction & control for site", subtitle = "lines based on\nlm(log(DPR+0.01)~mean_year*abs(Latitude)+(1|site)") + 
        theme(plot.title = element_text(size=9),
          plot.subtitle = element_text(size=8),
          legend.position="none")
      ggsave('Outputs/respEmail_Fig1C_N_interaction_siteControl_absLat_title.png', fC_title, , width = 3.2, height = 2.75) 
 
  # not used - plot predictions - controlled for population and using interaction hemisphere (not absolute)
    # predictions
      d = datapred[datapred$Belt %in% c('Arctic', 'North temperate'),]
      m <- lmer( log(DPR + 0.01) ~ mean_year*Belt + (1|pop_ID) , data = d) 
      nsim <- 5000
      bsim <- sim(m, n.sim=nsim) 
      v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
              apply(bsim@fixef, 2, quantile, prob=c(0.025,0.5,0.975))
      # values to predict for   
      newDa=data.frame(mean_year = seq(1968,2015, length.out=300), 
                      Belt = "Arctic")
      newDb=data.frame(mean_year = seq(1944,2013, length.out=300), 
                      Belt = "North temperate")
      newD = rbind(newDa,newDb)
      # exactly the model which was used has to be specified here 
        X <- model.matrix(~ mean_year*Belt,data=newD)  
                  
      # calculate predicted values and creditability intervals
        newD$pred <- exp(X%*%v)-0.01 # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(i in 1:nsim) predmatrix[,i] <- exp(X%*%bsim@fixef[i,])-0.01
            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
        arc=newD[newD$Belt == 'Arctic',]   
        NTemp=newD[newD$Belt == 'North temperate',]

        yearsA <- arc$mean_year
        yearsT <- NTemp$mean_year

        predarctic  <- arc$pred
        predtempN <- NTemp$pred
        predarcticUPR = arc$upr
        predarcticLWR = arc$lwr
        predtempNUPR = NTemp$upr
        predtempNLWR  = NTemp$lwr

      preda <- data.frame( years = yearsA, predarctic = predarctic,isarctic = TRUE )
      predaLwr <- data.frame( years = yearsA, predarctic = predarcticLWR,isarctic = TRUE )
      predaUpr <- data.frame( years = yearsA, predarctic = predarcticUPR,isarctic = TRUE )
      predo <- data.frame( years = yearsT, predtempN = predtempN, isarctic = FALSE )
      predoLwr <- data.frame( years = yearsT, predtempN = predtempNUPR, isarctic = FALSE )
      predoUpr <- data.frame( years = yearsT, predtempN = predtempNLWR, isarctic = FALSE )
    # plot
      figurebHem <- ggplot( datasummary, aes( x = decade*10+5, y = meandpr ) ) + aes( colour = isarctic) +
        labs(x = "Year", y = "Mean DPR") +
        ylim(0, 0.07)+ 
        geom_errorbar( aes(ymin = minusse, ymax = plusse ), width = 0.1, position = position_dodge(width = 0.25) )+
        geom_point( aes(size =n), position = position_dodge(width = 0.25)) +
        guides(size=guide_legend(title=expression(N[estimates])))+
        theme_bw() + 
        theme( panel.border = element_blank(), 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               #legend.position="none",
               axis.line = element_line(colour = 'black', size = 0.25),
               axis.ticks = element_line(colour = "black", size = 0.25),
               axis.text.x = element_text(size = 10),
               axis.text.y = element_text(size = 8),
               #axis.title.x=element_blank(),
               axis.title.y=element_text(size = 10) )
      figurebHem <- figurebHem + scale_color_manual(breaks = c( "FALSE", "TRUE"), values = c("darkolivegreen4", "deeppink" ), guide = FALSE)

      figurebHem <- figurebHem + 
          geom_line( aes(x = years, y = predarctic), preda  ) + 
          geom_line( aes(x = years, y = predarctic), predaLwr, lty = 3) + 
          geom_line( aes(x = years, y = predarctic), predaUpr, lty = 3) + 

          geom_line( aes(x = years, y = predtempN), predo  ) +
          geom_line( aes(x = years, y = predtempN), predoLwr, lty = 3) + 
          geom_line( aes(x = years, y = predtempN), predoUpr, lty = 3)
      ggsave('Outputs/respEmail_Fig1C_N_interaction_siteControl_ArcticVsNTemp.png',figurebHem, width = 3.5, height = 2.7)
# Figure R2 bottom -same as top but with TPR instead of DPR
  # prepare data summary 
    # means
       (dataRed <- datapred[c(which(datapred$Belt == "North temperate"), which(datapred$Belt == "Arctic")),])
       (datasummary <- dataRed %>% group_by( isarctic, decade) %>% 
        summarize( meandpr = ( mean((TPR)  )),
                   se = stderr( TPR ),
                   plusse =  ( mean((TPR)  ) + stderr( (TPR) )   ),
                   minusse =  ( mean((TPR)  ) -  stderr( (TPR) ) ),
                   sel  = stderr( (TPR) ),
                    n = length(DPR) ) )
  # plot Kubelka et al's oroginal Figure C but with TPR
    model.lm1 <- lm( TPR ~ poly(mean_year, 1) + Latitude  , data = datapred) 
    years <- seq(1944, 2016, by = 1)

    (newarctic <- data.frame( mean_year = years, Latitude = rep(69, 73) ))
    (newtempN  <- data.frame( mean_year = years, Latitude = rep(50.78, 73) ))

    predarctic  <- predict(model.lm1, newarctic)
    predtempN <- predict(model.lm1, newtempN) 

    preda <- data.frame( years = years, predarctic = predarctic, isarctic = TRUE )
    predo <- data.frame( years = years, predtempN = predtempN, isarctic = FALSE )

    figure <- ggplot( datasummary, aes( x = decade*10+5, y = meandpr ) ) + aes( colour = isarctic) +
      labs(x = "Year", y = "Total predation rate") +
      ylim(0, 0.73)+ 
      geom_errorbar( aes(ymin = minusse, ymax = plusse ), width = 0.1, position = position_dodge(width = 0.25) )+
      geom_point( position = position_dodge(width = 0.25)) +
      theme_bw() + 
      theme( panel.border = element_blank(), 
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             #legend.position="none",
             axis.line = element_line(colour = 'black', size = 0.25),
             axis.ticks = element_line(colour = "black", size = 0.25),
             axis.text.x = element_text(size = 10),
             axis.text.y = element_text(size = 8),
             #axis.title.x=element_blank(),
             axis.title.y=element_text(size = 10) )
    figure <- figure + scale_color_manual(breaks = c( "FALSE", "TRUE"), values = c("darkolivegreen4", "deeppink" ), guide = FALSE)

    figure <- figure + geom_line( aes(x = years, y = predarctic), preda  ) + geom_line( aes(x = years, y = predtempN), predo  )

    figure = figure + geom_text(x=1953, y=0.625, label="Arctic", size = 3, col ='grey30', adj = 0) + geom_point(x = 1950, y = 0.625,col = "deeppink")
    figure = figure + geom_text(x=1953, y=0.565, label="North temperate", col = "grey30", size = 3, adj = 0)+ geom_point(x = 1950, y = 0.565,col = "darkolivegreen4")
    
    #ggsave('Outputs/respEmail_Fig1C_original_TPR.png', figure, , width = 3.5, height = 2)
    
    f_title = figure + ggtitle(label ="Figure 1C as in Kubelka et al. 2019b, but on TPR", subtitle = "lines based on\nlm(TPR~poly(mean_year,1)+Latitude") + 
        theme(plot.title = element_text(size=9),
          plot.subtitle = element_text(size=8))
    ggsave('Outputs/Figure_R2_bottomLeft_Fig1C_original_TPR.png', f_title, , width = 3.2, height = 2.75) 
  # plot predictions - controlled for population and using interaction with linear latitude (not absolute)
    # prepare predictions
      m <- lmer( TPR ~ scale(mean_year)*scale(Latitude) + (1|pop_ID) , data = datapred)  
      summary(glht(m))
      m <- lmer( TPR ~ mean_year*Latitude + (1|pop_ID) , data = datapred) 
      nsim <- 5000
      bsim <- sim(m, n.sim=nsim) 
      v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
              apply(bsim@fixef, 2, quantile, prob=c(0.025,0.5,0.975))
          # values to predict for   
          newD=data.frame(mean_year = seq(1968,2015, length.out=300),
                          Latitude = 69)
           
          # exactly the model which was used has to be specified here 
            X <- model.matrix(~ mean_year*Latitude,data=newD)  
                      
          # calculate predicted values and creditability intervals
            newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
                newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
            arc=newD   

          # values to predict for   
          newD=data.frame(mean_year = seq(1944,2013, length.out=300),
                          Latitude = 50.78)
           
          # exactly the model which was used has to be specified here 
            X <- model.matrix(~ mean_year*Latitude,data=newD)  
                      
          # calculate predicted values and creditability intervals
            newD$pred <- X%*%v
            predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
            for(i in 1:nsim) predmatrix[,i] <- (X%*%bsim@fixef[i,])
                newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
                newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
            NTemp=newD  

        yearsA <- arc$mean_year
        yearsT <- NTemp$mean_year

        predarctic  <- arc$pred
        predtempN <- NTemp$pred
        predarcticUPR = arc$upr
        predarcticLWR = arc$lwr
        predtempNUPR = NTemp$upr
        predtempNLWR  = NTemp$lwr

      preda <- data.frame( years = yearsA, predarctic = predarctic,isarctic = TRUE )
      predaLwr <- data.frame( years = yearsA, predarctic = predarcticLWR,isarctic = TRUE )
      predaUpr <- data.frame( years = yearsA, predarctic = predarcticUPR,isarctic = TRUE )
      predo <- data.frame( years = yearsT, predtempN = predtempN, isarctic = FALSE )
      predoLwr <- data.frame( years = yearsT, predtempN = predtempNUPR, isarctic = FALSE )
      predoUpr <- data.frame( years = yearsT, predtempN = predtempNLWR, isarctic = FALSE )
    # plot
      figurebLinLat <- ggplot( datasummary, aes( x = decade*10+5, y = meandpr ) ) + aes( colour = isarctic) +
        labs(x = "Year", y = "Total predation rate") +
        ylim(0, 0.73)+ 
        geom_errorbar( aes(ymin = minusse, ymax = plusse ), width = 0.1, position = position_dodge(width = 0.25) )+
        geom_point( aes(size =n), position = position_dodge(width = 0.25)) +
        guides(size=guide_legend(title=expression(N[estimates])))+
        theme_bw() + 
        theme( panel.border = element_blank(), 
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(), 
               #legend.position="none",
               axis.line = element_line(colour = 'black', size = 0.25),
               axis.ticks = element_line(colour = "black", size = 0.25),
               axis.text.x = element_text(size = 10),
               axis.text.y = element_text(size = 8),
               #axis.title.x=element_blank(),
               axis.title.y=element_text(size = 10) )
      figurebLinLat <- figurebLinLat + scale_color_manual(breaks = c( "FALSE", "TRUE"), values = c("darkolivegreen4", "deeppink" ), guide = FALSE)
      figurebLinLat <- figurebLinLat + 
        geom_line( aes(x = years, y = predarctic), preda  ) + 
        geom_line( aes(x = years, y = predarctic), predaLwr, lty = 3) + 
        geom_line( aes(x = years, y = predarctic), predaUpr, lty = 3) + 

        geom_line( aes(x = years, y = predtempN), predo  ) +
        geom_line( aes(x = years, y = predtempN), predoLwr, lty = 3) + 
        geom_line( aes(x = years, y = predtempN), predoUpr, lty = 3)
      #ggsave('Outputs/respEmail_Fig1C_N_interaction_siteControl_latNotabs_TPR.png',figurebLinLat, width = 3.2, height = 2.75)

      fC_title = figurebLinLat + ggtitle(label ="Test for interaction & control for site", subtitle = "lines based on\nlm(TPR~mean_year*Latitude+(1|site)") + 
        theme(plot.title = element_text(size=9),
          plot.subtitle = element_text(size=8),
          legend.position="none")
      ggsave('Outputs/Figure_R2_bottomLeft_Fig1C_interaction_siteControl_latNotabs_PR.png', fC_title, , width = 3.2, height = 2.75) 
              
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
  #  [1] multcomp_1.4-13  TH.data_1.0-10   survival_3.1-12  mvtnorm_1.1-1    dplyr_1.0.1     
  #  [6] gridExtra_2.3    ggplot2_3.3.2    AICcmodavg_2.3-1 lmerTest_3.1-3   MuMIn_1.43.17   
  # [11] arm_1.11-2       lme4_1.1-23      Matrix_1.2-18    MASS_7.3-51.6   
  # 
  # loaded via a namespace (and not attached):
  #  [1] VGAM_1.1-3          splines_4.0.2       Formula_1.2-3       statmod_1.4.34     
  #  [5] sp_1.4-2            stats4_4.0.2        latticeExtra_0.6-29 numDeriv_2016.8-1.1
  #  [9] pillar_1.4.6        backports_1.1.8     lattice_0.20-41     glue_1.4.2         
  # [13] digest_0.6.25       RColorBrewer_1.1-2  checkmate_2.0.0     minqa_1.2.4        
  # [17] colorspace_1.4-1    sandwich_2.5-1      htmltools_0.5.0     plyr_1.8.6         
  # [21] pkgconfig_2.0.3     raster_3.3-13       purrr_0.3.4         xtable_1.8-4       
  # [25] scales_1.1.1        jpeg_0.1-8.1        htmlTable_2.0.1     tibble_3.0.3       
  # [29] farver_2.0.3        generics_0.0.2      ellipsis_0.3.1      withr_2.2.0        
  # [33] nnet_7.3-14         magrittr_1.5        crayon_1.3.4        nlme_3.1-148       
  # [37] foreign_0.8-80      tools_4.0.2         data.table_1.13.0   lifecycle_0.2.0    
  # [41] stringr_1.4.0       munsell_0.5.0       cluster_2.1.0       compiler_4.0.2     
  # [45] rlang_0.4.7         grid_4.0.2          nloptr_1.2.2.2      rstudioapi_0.11    
  # [49] htmlwidgets_1.5.1   base64enc_0.1-3     labeling_0.3        boot_1.3-25        
  # [53] gtable_0.3.0        codetools_0.2-16    abind_1.4-5         R6_2.4.1           
  # [57] zoo_1.8-8           knitr_1.29          Hmisc_4.4-0         stringi_1.5.3      
  # [61] parallel_4.0.2      unmarked_1.0.1      Rcpp_1.0.5          vctrs_0.3.2        
  # [65] rpart_4.1-15        acepack_1.4.1       png_0.1-7           tidyselect_1.1.0   
  # [69] xfun_0.16           coda_0.19-3        