# TOOLS and DATA manipulation from Kubelka et al
  rm( list = ls() )
  # load packages
    library( ape )
    library( coxme )
    library( phytools )
    library(ggplot2)
    library(mgcv)
    require(arm)
    require(here)
    require(multcomp)
    require(ggplot2)
    require(effects)
    require(plyr)
    require(glue)
    require(RColorBrewer)
    library(lme4)
    library(lmerTest)
  # FUNCTIONS
    # function for calculating distances -
    geodetic <- function(l1, t1, l2, t2) {
      
      l1 <- l1 / 360.0 * 2.0 * pi
      l2 <- l2 / 360.0 * 2.0 * pi
      t1 <- t1 / 360.0 * 2.0 * pi
      t2 <- t2 / 360.0 * 2.0 * pi
      
      dist <- 6371.0 * acos( sin(t1) * sin(t2) + cos(t1) * cos(t2) * cos(l2 - l1) )
      
      return(dist)
      
    }

    # function for a distance matrix
    dist.mat <- function(lat, lon, rwnms) {
      n <- length( lat )
      mat <- matrix(0, n, n )
      
      for(i in 1:n) {
        for( j in 1:n) {
          mat[i, j] <- geodetic( lon[i], lat[i], lon[j], lat[j] )
          if( is.na( mat[i,j]) == TRUE ) mat[i,j] <- 0    # Nasty
        }
      }
      
      mdist <- geodetic(-90,-90, 90,90)
      mat <- mdist - mat
      diag(mat) <- mdist
      mat <- mat / mdist
      rownames(mat) <- rwnms
      return(mat)
      
    }

    # Add in the missing species
    bind.tip<-function(tree,tip.label,edge.length=NULL,where=NULL){
      if(is.null(where)) where<-length(tree$tip)+1
      tip<-list(edge=matrix(c(2,1),1,2),
                tip.label=tip.label,
                edge.length=edge.length,
                Nnode=1)
      class(tip)<-"phylo"
      obj<-bind.tree(tree,tip,where=where)
      return(obj)
    }

    addInTip <- function( phylo, where, newname) {
      idx <- which( phylo$tip == where )
      np <- nodepath( phylo)[[idx]]
      to <- np[ length(np) -1 ]
      ed <- which( tree$edge[,2] == idx )
      newphylo <- bind.tip( phylo, newname, edge.length = phylo$edge.length[ed], where = to )
      return( newphylo)
    }
  # load and prepare trees  
    trees <- read.tree( "Data/trees2.phy" )
    tree <- trees[[42]]  # just choose one to start with
    idx <- which( tree$tip == "Charadrius_alexandrinus")
    np <- nodepath(tree)[[idx]]
    to <- np[ length( np ) - 1] 
    ed <- which( tree$edge[,2] == idx)

    tree1 <- bind.tip( tree, "Charadrius_nivosus", edge.length = tree$edge.length[ed], where = to )
    plotTree(tree1, node.numbers = TRUE, fsize = 0.5)
    tree2 <- addInTip( tree1,"Gallinago_gallinago", "Gallinago_delicata" )
  # load and prepare data
    data <- read.csv( "Data/DATApopulations.csv", h = T, sep=";")	
	  datapred <- read.csv("Kubelka et al. 2019 Response_Data, R codes and supporting information/DATApopulations.csv", sep =";", h = T)

    # prepare subsets
      datapred$pop_ID <-paste(datapred$Longitude, datapred$Latitude)
      datapred$hemis <-  datapred$Latitude > 0 
      datapred$arctic <- datapred$Belt == "Arctic"
      sub97 <- datapred[which(datapred$DPR_trans == "NO"),]
      tab <- table( sub97$pop_ID )
      idx <- which(tab > 1)
      inclPops <- names(tab)[idx]
      retain <- which( sub97$pop_ID %in% inclPops == TRUE)
      subsub97 <- sub97[retain,]

    # Create latitude and site variables
    	data$lat_abs = abs(data$Latitude)
    	data$NS <- data$Latitude > 0
    	data$site = paste(data$Latitude,data$Longitude)
  	
    #  Create phylogenetic geographic and PI matrices
      phyloMat <- vcv.phylo( tree2 )
      phyloMat <- phyloMat / max( phyloMat )
      distanceMatrix <- dist.mat( data$Latitude, data$Longitude, data$species ) 
      diag( distanceMatrix ) <- diag(distanceMatrix) + 0.01
      distanceMatrix <- distanceMatrix / 1.01

      idxMat <- apply( phyloMat, 1:2, function(i) as.numeric( i > 0) )


      I0 <- diag(1, dim(data)[1] )
      rownames(I0) <- colnames(I0) <- data$species 
       I <- diag(1 / data$N_nests )
       
      I <- I / max(I)
      rownames(I) <- colnames(I) <- data$species 

    # Form a matrix representing species identity
      m <- as.matrix( model.matrix( ~ species-1, data = data) )
      S <- crossprod(t(m), t(m) )
      rownames(S) <- colnames(S) <- data$species 

# Note 2 - scrutinizing Kubelka et al's response models about global increase in predation rates over time
  # M2 - model used in Kubelka et al. (2019a) as a replication of our model (lmer( log(DPR + 0.01) ~ log(N_nests) + Latitude + mean_year + (1|pop_ID))
    m2 <- lmer( log(DPR + 0.01) ~ mean_year + (1|pop_ID) , data = subsub97) 
    summary(glht(m2))
  
  # Picture A - M2 fitted to all directly calculated estimates (as we did in our Comment)
    m2c_all <- lmer( log(DPR + 0.01) ~  mean_year  + (1|pop_ID) , data = sub97) 
    summary(m2c_all)
    summary(glht(m2c_all))
  
  # Picture B - M2 controlled for number of nests to generate an estimate and for latitude of the population
    m2c <- lmer( log(DPR + 0.01) ~ log(N_nests) + Latitude + mean_year + (1|pop_ID) , data = subsub97) 
    summary(m2c)
    summary(glht(m2c))

  # Picture C - M2 controlled for number of nests to generate an estimate and including interaction with latitude as in our original Figure 1F (referred to by Kubelka et al. 2019) 
    m2c_asF3E <- lmer( log(DPR + 0.01) ~ log(N_nests) + scale(mean_year)*poly(Latitude,3) + (1|pop_ID) , data = subsub97) 
    summary(m2c_asF3E)
    summary(glht(m2c_asF3E))  

    m2c_asF3E <- lmer( log(DPR + 0.01) ~ log(N_nests) + scale(mean_year)*scale(Latitude) + (1|pop_ID) , data = subsub97) 
    summary(m2c_asF3E)
    summary(glht(m2c_asF3E))

# Note 5 - Fig DPR original vs log
  source('R/Constants_Functions.R')
  d = data
  d$Belt = as.factor(d$Belt) # define site
  d$site = paste(d$Latitude,d$Longitude) # define site
  d$lat_abs = abs(d$Latitude) # abs latitude
  d$ln_N_nests = log(d$N_nests)
  d$hemisphere =as.factor(ifelse(d$Latitude > 0, "Northern", "Southern"))
  d$genus = gsub("\\_.*","",d$species)
  d$DPR_trans[which(d$DPR_trans == 'YES' & d$source_id=="209")] = "NO"# source_id 209 Schekkerman et al. 1998 (Calidris 

  # transformed
    dd =d[d$DPR_trans == 'NO' & d$Belt %in% c('Arctic','North temperate'),]
    nrow(dd)
    #summary(dd$Belt) 
    m = lmer(log(DPR) ~ ln_N_nests + poly(mean_year,2)*Belt+(1|site)+(1|species),  data =dd )
      nsim <- 5000
      bsim <- sim(m, n.sim=nsim)  
      
      # coefficients
      v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          apply(bsim@fixef, 2, quantile, prob=c(0.025,0.5,0.975))
      # values to predict for   
      newD=data.frame(ln_N_nests = mean(dd$ln_N_nests),mean_year = seq(min(dd$mean_year),max(dd$mean_year), length.out=900), Belt = c('Arctic', 'North temperate'))
      
            
      # exactly the model which was used has to be specified here 
        X <- model.matrix(~ ln_N_nests + poly(mean_year,2)*Belt,data=newD)  
                  
      # calculate predicted values and creditability intervals
        newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
        newD$line_col = col_$line_col[match(newD$Belt, col_$Belt)]
    dprC=newD   
  
  # original
    dd =d[d$DPR_trans == 'NO' & d$Belt %in% c('Arctic','North temperate'),]
    nrow(dd)
    #summary(dd$Belt) 
    m = lmer(DPR ~ ln_N_nests + poly(mean_year,2)*Belt+(1|site)+(1|species),  data =dd )
      nsim <- 5000
      bsim <- sim(m, n.sim=nsim)  
      
      # coefficients
      v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
          apply(bsim@fixef, 2, quantile, prob=c(0.025,0.5,0.975))
      # values to predict for   
      newD=data.frame(ln_N_nests = mean(dd$ln_N_nests),mean_year = seq(min(dd$mean_year),max(dd$mean_year), length.out=900), Belt = c('Arctic', 'North temperate'))
      
            
      # exactly the model which was used has to be specified here 
        X <- model.matrix(~ ln_N_nests + poly(mean_year,2)*Belt,data=newD)  
                  
      # calculate predicted values and creditability intervals
        newD$pred <- X%*%v # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(i in 1:nsim) predmatrix[,i] <- X%*%bsim@fixef[i,]
            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
        newD$line_col = col_$line_col[match(newD$Belt, col_$Belt)]
    dprCorig=newD

  # PLOT horizontal
    if(PNG == TRUE){
    #Cairo(file=glue(here::here(),"/Outputs/original_vs_ln-transformed_horizontalTEST.png"), type="png", units="in",  width=1.6*2,height=1.8, dpi=300)
    png(glue(here::here(),"/Outputs/original_vs_log-transformed_horizontal_base.png"), width=1.6*2,height=1.8,units="in",res=600) 
     }else{dev.new(width=1.6*2,height=1.8)}
    
    par(oma = c(1.0, 1.5, 0.55, 0.25),mgp=c(1.2,0.15,0), 
        ps=12,font.main = 1, las=1, lwd = 0.5, tcl=-0.05,
        cex=1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5,
        col.lab="black", col.main="162", fg="grey45",
        bty="n",xpd=TRUE) #
  
    layout(mat = matrix(c(1,2),ncol = 2))
     

    pp = dprCorig

    par(mar=c(0.1,0.5,0,0.1),ps=12, cex=1, font.main = 1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.05,bty="n",xpd=FALSE)
   
    plot(pp$pred~pp$mean_year, pch=19,xlim=c(1944,2016), ylim=c(0,0.07),ylab=NA,xlab=NA, xaxt='n',yaxt='n', type='n')#xaxs="i",yaxs="i")
    
    axis(1, at=c(1944,1960,1980,2000,2016),labels=c(1944,1960,1980,2000,2016),cex.axis=0.5,mgp=c(0,-0.35,0), lwd = 0.5)
    axis(2, at=seq(0,0.07, by = 0.01), labels=c('0','0.01','0.02','0.03','0.04','0.05','0.06','0.07'), lwd = 0.5)
    mtext('Daily nest predation',side=2,line=1.1, cex=0.6, las=3, col = 'black')
    #text(x = 2016+5, y =0.07*1.01, labels= expression(bold("Based on model with predation on")), col='black', cex = 0.6, xpd = TRUE)
    mtext('Year',side=1,line=0.2, cex=0.6, las=1, col = 'black', at = 2016+10 )
    mtext(expression(bold("Based on model with predation rate on:")), side = 3, line = -0.2, col='black', cex = 0.6, xpd = TRUE, at = 2016+10)
    mtext("original scale", side = 3, line = -0.8, col='black', cex = 0.6, xpd = TRUE, at =mean(c(1944,2016)))
    #text(x = 2016-30, y =0.07*0.98, labels= expression(bold("original scale")), col='black', cex = 0.6, xpd = TRUE)
    for(i in unique(pp$Belt)){
    #i ="Arctic"
    print(i)
    ppi = pp[pp$Belt==i,]
    polygon(c(ppi$mean_year, rev(ppi$mean_year)), c((ppi$lwr)-0.01, 
      rev((ppi$upr)-0.01)), border=NA, col = adjustcolor(ppi$line_col,alpha.f = 0.1))#adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
    lines(ppi$mean_year, (ppi$pred)-0.01, col=ppi$line_col,lwd=1)
            } 
    
    pp = dprC
    #par(mar=c(0.1,0.5,0,0.1),ps=12, cex=1, font.main = 1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.05,bty="n",xpd=FALSE)
    par(mar=c(0.1,0.4,0,0.2),ps=12, cex=1, font.main = 1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.05,bty="n",xpd=FALSE)
   
    plot(exp(pp$pred)~pp$mean_year, pch=19,xlim=c(1944,2016), ylim=c(0,0.07),ylab=NA,xlab=NA, xaxt='n',yaxt='n', type='n')#xaxs="i",yaxs="i")
    
    axis(1, at=c(1944,1960,1980,2000,2016),labels=c(1944,1960,1980,2000,2016),cex.axis=0.5,mgp=c(0,-0.35,0), lwd = 0.5)
    #mtext('Year',side=1,line=0.2, cex=0.6, las=1, col = 'black')
    axis(2, at=seq(0,0.07, by = 0.01), labels=NA, lwd = 0.5)
    #mtext('Daily nest predation',side=2,line=1.1, cex=0.6, las=3, col = 'black')
    mtext("log-scale", side = 3, line = -0.8, col='black', cex = 0.6, xpd = TRUE, at =mean(c(1944,2016)))
    
    #text(x = 2016-30, y =0.07*0.98, labels= expression(bold("ln-transformed")), col='black', cex = 0.7,  xpd = TRUE)
    for(i in unique(pp$Belt)){
      #i ="Arctic"
      print(i)
      ppi = pp[pp$Belt==i,]
      polygon(c(ppi$mean_year, rev(ppi$mean_year)), c(exp(ppi$lwr)-0.01, 
        rev(exp(ppi$upr)-0.01)), border=NA, col = adjustcolor(ppi$line_col,alpha.f = 0.1))#adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
      lines(ppi$mean_year, exp(ppi$pred)-0.01, col=ppi$line_col,lwd=1)
              } 
    if(PNG == TRUE){dev.off()}
  
  quartz.save(glue(here::here(),"/Outputs/original_vs_ln-transformed_horizontal_quartz.png"), type = "png") 
  
 # PLOT vertical
    dev.new(width=1.6,height=2*1.45)

    par(oma = c(1.0, 1.5, 0.55, 0.25),mgp=c(1.2,0.15,0), 
        ps=12,font.main = 1, las=1, lwd = 0.5, tcl=-0.05,
        cex=1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5,
        col.lab="black", col.main="162", fg="grey45",
        bty="n",xpd=TRUE) #
  
    layout(mat = matrix(c(1,2),nrow = 2))
     

    pp = dprCorig
    par(mar=c(0.1,0.5,0,0.1),ps=12, cex=1, font.main = 1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.05,bty="n",xpd=FALSE)
   
    plot(pp$pred~pp$mean_year, pch=19,xlim=c(1944,2016), ylim=c(0,0.07),ylab=NA,xlab=NA, xaxt='n',yaxt='n', type='n')#xaxs="i",yaxs="i")
    
    axis(1, at=c(1944,1960,1980,2000,2016),labels=NA,cex.axis=0.5,mgp=c(0,-0.35,0), lwd = 0.5)
    #mtext('Year',side=1,line=0.2, cex=0.6, las=1, col = 'black')
    axis(2, at=seq(0,0.07, by = 0.01), labels=c('0','0.01','0.02','0.03','0.04','0.05','0.06','0.07'), lwd = 0.5)
    mtext('Daily nest predation',side=2,line=1.1, cex=0.6, las=3, col = 'black')
    text(x = 2016-30, y =0.07*0.98, labels= expression(bold("original scale")), col='black', cex = 0.7, xpd = TRUE)
    for(i in unique(pp$Belt)){
    #i ="Arctic"
    print(i)
    ppi = pp[pp$Belt==i,]
    polygon(c(ppi$mean_year, rev(ppi$mean_year)), c((ppi$lwr)-0.01, 
      rev((ppi$upr)-0.01)), border=NA, col = adjustcolor(ppi$line_col,alpha.f = 0.1))#adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
    lines(ppi$mean_year, (ppi$pred)-0.01, col=ppi$line_col,lwd=1)
            } 
    
    pp = dprC
    par(mar=c(0.1,0.5,0,0.1),ps=12, cex=1, font.main = 1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.05,bty="n",xpd=FALSE)
     
    plot(exp(pp$pred)~pp$mean_year, pch=19,xlim=c(1944,2016), ylim=c(0,0.07),ylab=NA,xlab=NA, xaxt='n',yaxt='n', type='n')#xaxs="i",yaxs="i")
    
    axis(1, at=c(1944,1960,1980,2000,2016),labels=c(1944,1960,1980,2000,2016),cex.axis=0.5,mgp=c(0,-0.35,0), lwd = 0.5)
    mtext('Year',side=1,line=0.2, cex=0.6, las=1, col = 'black')
    axis(2, at=seq(0,0.07, by = 0.01), labels=c('0','0.01','0.02','0.03','0.04','0.05','0.06','0.07'), lwd = 0.5)
    #mtext('Daily nest predation',side=2,line=1.1, cex=0.6, las=3, col = 'black')
    text(x = 2016-30, y =0.07*0.98, labels= expression(bold("ln-transformed")), col='black', cex = 0.7,  xpd = TRUE)
    for(i in unique(pp$Belt)){
      #i ="Arctic"
      print(i)
      ppi = pp[pp$Belt==i,]
      polygon(c(ppi$mean_year, rev(ppi$mean_year)), c(exp(ppi$lwr)-0.01, 
        rev(exp(ppi$upr)-0.01)), border=NA, col = adjustcolor(ppi$line_col,alpha.f = 0.1))#adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
      lines(ppi$mean_year, exp(ppi$pred)-0.01, col=ppi$line_col,lwd=1)
              } 
            
 quartz.save(glue(here::here(),"/Outputs/original_vs_ln-transformed_horizontal.png"), type = "png")

# Note 7 - scrutinizing Kubelka et al's response models about global increase in predation rates over time
  # M3 not controlled for multiple sampling per site
      M3 <- lm( log( DPR + 0.01) ~ log(N_nests) + mean_year +  Latitude , data = sub97) 
      summary(M3)
  # M4 adds the controls for multiple sampling per site
      M4 <- lmer( log(DPR + 0.01) ~ log(N_nests) + mean_year  +  Latitude + (1|pop_ID) , data = sub97) 
      summary(M4)
        m <- lmer( log(DPR + 0.01) ~ log(N_nests) + mean_year  + (1|pop_ID) , data = sub97) # simplified

  # M4 not controlled for number of nests for a given estimate and latitude 
      M4<- lmer( log(DPR + 0.01) ~ mean_year  + (1|pop_ID) , data = sub97) 
      summary(M4)


  tab <- table( sub97$pop_ID )
  idx <- which(tab > 1)
  inclPops <- names(tab)[idx]


  retain <- which( sub97$pop_ID %in% inclPops == TRUE)
  subsub97 <- sub97[retain,]
  model.lmer <- lmer( log(DPR + 0.01) ~ mean_year + (1|pop_ID) , data = subsub97) 
  summary(glht(model.lmer))

  m <- lmer( log(DPR + 0.01) ~ log(N_nests) + Latitude + mean_year + (1|pop_ID) , data = subsub97) 
  summary(glht(m))
  m <- lmer( log(DPR + 0.01) ~ log(N_nests) + mean_year  + (1|pop_ID) , data = subsub97) 
  m <- lmer( log(DPR + 0.01) ~ Latitude + mean_year  + (1|pop_ID) , data = subsub97) 
# Note 9 - results based on the mean year estimates from less than 10 years of data
  d = data
  dd =d[d$DPR_trans == 'NO' & d$Belt %in% c('Arctic','North temperate') & data$years_nr<10,]
  nrow(dd)
  m = lmer( log(DPR + 0.01) ~ log(N_nests) + mean_year  + Latitude + (1|site) , data = dd)
  m = lmer( log(DPR + 0.01) ~  mean_year  + (1|site) , data = dd)
  summary(m)
  summary(glht(m))