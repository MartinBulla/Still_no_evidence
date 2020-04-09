rm( list = ls() )
#setwd('/Users/martinbulla/Dropbox/Science/ms_published/Kubelka_et_al_rebuttal/dryad/')
#setwd('C:/Users/mbulla/Documents/Dropbox/Science/MS/Kubelka_et_al_rebuttal/dryad/')

#-----------------------
# TOOLS and DATA manipulation from Kubelka et al
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

# distribution of sites over globe and sample sizes 
  ggplot(data, aes(y = Latitude, x = Longitude)) + geom_jitter(aes(size = N_nests), alpha = 0.5)
  quartz.save(glue(here::here(),"/Outputs/N-nest_per_estimate_geography.png"), type = "png")
  ggsave(file = glue(here::here(),"/Outputs/N-nest_per_estimate_geography.png", dpi = 300, width = 10, height = 8, units = "cm") 
  
  xx = ddply(data, .(Latitude, Longitude), summarise, N_estimates = length(Latitude))
  ggplot(xx, aes(y = Latitude, x = Longitude)) + geom_point(aes(size = N_estimates), alpha = 0.5) + scale_size(limits = c(1, 15)) + scale_size(name = "# of estimates")
  quartz.save(glue(here::here(),"/Outputs/N-estimates_per_site_geography.png"), type = "png")
  ggsave(file=glue(here::here(),"/Outputs/N-estimates_per_site_geography.png"), dpi = 300, width = 10, height = 8, units = "cm")

  
  xx = ddply(data, .(Latitude, Longitude, period_orig), summarise, N_estimates = length(Latitude))
  ggplot(xx, aes(y = Latitude, x = Longitude, col = period_orig)) + 
      geom_point(aes(size = N_estimates), alpha = 0.5) + 
      scale_size(name = "# of estimates", breaks = c(1, 5, 10), limits = c(1, 15))+ 
      scale_color_manual(values = c('red','blue'), name = "Mean year of the study", labels =c("\u2265 2000", "< 2000"))
  ggsave(file=glue(here::here(),"/Outputs/N-estimates_per_site_geography_period.png"), dpi = 300, width = 12, height = 8, units = "cm")

  ggplot(data, aes(x = abs(Latitude), y = TPR, col = mean_year)) + geom_jitter(aes(size = N_nests), alpha = 0.5)
  ggplot(data, aes(x = abs(Latitude), y = DPR, col = mean_year)) + geom_jitter(aes(size = N_nests), alpha = 0.5)
    

  xy = ddply(data, .(Latitude, Longitude, mean_year, period_orig, Belt), summarise, N_estimates = length(Latitude), med_TPR = median(TPR), avg_TPR = mean(TPR), med_DPR = median(DPR), avg_DPR = mean(DPR))
 
  ggplot(xy, aes(y = abs(Latitude), x = med_TPR, col = mean_year)) + geom_point(aes(size = N_estimates), alpha = 0.8) + scale_size(limits = c(1, 15)) + scale_colour_gradientn(colours = rev(brewer.pal(11, "Spectral")))
    ggsave(file=glue(here::here(),"/Outputs/N-estimates_per_site_geography_medTPR_year.png"), dpi = 300, width = 10, height = 8, units = "cm")

  ggplot(xy, aes(y = abs(Latitude), x = avg_TPR, col = mean_year)) + geom_point(aes(size = N_estimates), alpha = 0.8) + scale_size(limits = c(1, 15)) + scale_colour_gradientn(colours = rev(brewer.pal(11, "Spectral")))
    ggsave(file=glue(here::here(),"/Outputs/N-estimates_per_site_geography_avgTPR_year.png"), dpi = 300, width = 10, height = 8, units = "cm")

  ggplot(xy, aes(y = abs(Latitude), x = med_TPR, col = period_orig)) + geom_point(aes(size = N_estimates), alpha = 0.5) + scale_size(limits = c(1, 15)) + scale_color_manual(values = c('red','blue'))
    ggsave(file=glue(here::here(),"/Outputs/N-estimates_per_site_geography_medTPR_period.png"), dpi = 300, width = 10, height = 8, units = "cm")

  ggplot(xy, aes(y = abs(Latitude), x = med_DPR, col = mean_year)) + geom_point(aes(size = N_estimates), alpha = 0.8) + scale_size(limits = c(1, 15)) + scale_colour_gradientn(colours = rev(brewer.pal(11, "Spectral")))
    ggsave(file=glue(here::here(),"/Outputs/N-estimates_per_site_geography_medDPR_year.png"), dpi = 300, width = 10, height = 8, units = "cm")

  ggplot(xy, aes(y = abs(Latitude), x = avg_DPR, col = mean_year)) + geom_point(aes(size = N_estimates), alpha = 0.8) + scale_size(limits = c(1, 15)) + scale_colour_gradientn(colours = rev(brewer.pal(11, "Spectral")))
    ggsave(file=glue(here::here(),"/Outputs/N-estimates_per_site_geography_avgDPR_year.png"), dpi = 300, width = 10, height = 8, units = "cm")

  ggplot(xy, aes(y = abs(Latitude), x = med_DPR, col = period_orig)) + geom_point(aes(size = N_estimates), alpha = 0.5) + scale_size(limits = c(1, 15)) + scale_color_manual(values = c('red','blue'))
    ggsave(file=glue(here::here(),"/Outputs/N-estimates_per_site_geography_medDPR_period.png"), dpi = 300, width = 10, height = 8, units = "cm")

# DATA summaries
  xx = data.frame(table(data$site))
  table(xx$Freq)

  xa = data.frame(table(data$site[data$Belt == 'Arctic']))
  table(xa$Freq)

  summary(factor(data$years_nr)) 
  nrow(data[data$years_nr>9,])

  ggplot(data, aes(x = Belt, y = years_nr)) +
    geom_violin(trim = FALSE, fill = '#A4A4A4') + 
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=1) +  
    geom_boxplot(width=0.1)+
    theme_minimal()

  ggplot(data, aes(x = Belt, y = years_nr, fill = period_orig)) +
    #geom_violin(trim = FALSE) + 
    geom_boxplot()+
    theme_minimal()  

  ggplot(data, aes(x = years_nr, y = TPR, col = period_orig)) +
    geom_point()+ stat_smooth()+
    theme_minimal()    
 
 ggplot(data, aes(x = mean_year, y = years_nr)) +
    geom_point()+ stat_smooth()+
    theme_minimal() 

# CHECK INDEPENDENCE OF RESIDUALS from Kubelka et al's original Table S2 models
  # TABLE S2A model - General
    m0 <- lmekin( log(DPR) ~ (1|species ) + mean_year + log( N_nests) , varlist = list( I, phyloMat , distanceMatrix ), data = data )
    #m0 <- lmekin( log(DPR) ~ (1|species ) + Belt+mean_year + log( N_nests) , varlist = list( I, phyloMat , distanceMatrix ), data = data )

    # We show the lack of control for pseudoreplication due to multiple data points per site, and especially for sites with multiple data
    data$res = resid(m0)
    mr <- lmer(res ~ 1+(1|site),  data = data)
    summary(mr) # site explain 56% of var in residuals
       l=data.frame(summary(mr)$varcor)
    	 l = l[is.na(l$var2),]
    	 l$var1 = ifelse(is.na(l$var1),"",l$var1)
    	 l$pred = paste(l$grp,l$var1)
    	ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
    	ri$estimate_r = paste(ri$estimate_r,'%',sep='')
    	ri

  # TABLE S2B model - South Temperate
    data <- read.csv( "South_temperate.csv", sep=";")
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

    # Create a latitude variable

    data$NS <- data$Latitude > 0


    m0 <- lmekin( log(DPR) ~ (1|species ) + mean_year + log( N_nests) , varlist = list( I, phyloMat , distanceMatrix ), data = data )
    
    # We show the lack of control for pseudoreplication due to multiple data points per site, and especially for sites with multiple data
    data$res = resid(m0)
    mr <- lmer(res ~ 1+(1|site),  data = data)
    summary(mr) # site explain 56% of var in residuals
       l=data.frame(summary(mr)$varcor)
       l = l[is.na(l$var2),]
       l$var1 = ifelse(is.na(l$var1),"",l$var1)
       l$pred = paste(l$grp,l$var1)
      ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
      ri$estimate_r = paste(ri$estimate_r,'%',sep='')
      ri

  # TABLE S2C model - South_tropics area (-30-0° latitude)

    data <- read.csv( "South_tropics.csv", sep=";")
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

    # Create a latitude variable

    data$NS <- data$Latitude > 0


    m0 <- lmekin( log(DPR) ~ (1|species ) + mean_year + log( N_nests) , varlist = list( I, phyloMat , distanceMatrix ), data = data )
    data$res = resid(m0)
    mr <- lmer(res ~ 1+(1|site),  data = data)
    summary(mr) # site explain 56% of var in residuals
       l=data.frame(summary(mr)$varcor)
       l = l[is.na(l$var2),]
       l$var1 = ifelse(is.na(l$var1),"",l$var1)
       l$pred = paste(l$grp,l$var1)
      ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
      ri$estimate_r = paste(ri$estimate_r,'%',sep='')
      ri

  # TABLE S2D North_tropics area (0-30° latitude)
    data <- read.csv( "North_tropics.csv", sep=";")
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

    # Create a latitude variable

    data$NS <- data$Latitude > 0


   m0 <- lmekin( log(DPR) ~ (1|species ) + mean_year + log( N_nests) , varlist = list( I, phyloMat , distanceMatrix ), data = data )
    data$res = resid(m0)
      mr <- lmer(res ~ 1+(1|site),  data = data)
      summary(mr) # site explain 56% of var in residuals
         l=data.frame(summary(mr)$varcor)
         l = l[is.na(l$var2),]
         l$var1 = ifelse(is.na(l$var1),"",l$var1)
         l$pred = paste(l$grp,l$var1)
        ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
        ri$estimate_r = paste(ri$estimate_r,'%',sep='')
        ri

  # TABLE S2E North_temperate area (0-30° latitude)

    data <- read.csv( "North_temperate.csv", sep=";")
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

    # Create a latitude variable

    data$NS <- data$Latitude > 0


    m0 <- lmekin( log(DPR) ~ (1|species ) + mean_year + log( N_nests) , varlist = list( I, phyloMat , distanceMatrix ), data = data )
    data$res = resid(m0)
      mr <- lmer(res ~ 1+(1|site),  data = data)
      summary(mr) # site explain 56% of var in residuals
         l=data.frame(summary(mr)$varcor)
         l = l[is.na(l$var2),]
         l$var1 = ifelse(is.na(l$var1),"",l$var1)
         l$pred = paste(l$grp,l$var1)
        ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
        ri$estimate_r = paste(ri$estimate_r,'%',sep='')
        ri

  # TABLE S2F Arctic (60-78° latitude)

    data <- read.csv( "Arctic.csv", sep=";")
    data$site = paste(data$Latitude,data$Longitude)
    # Create phylogenetic geographic and PI matrices
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


    m0 <- lmekin( log(DPR) ~ (1|species) + mean_year + log( N_nests), varlist = list( I, phyloMat, distanceMatrix), data = data )
      data$res = resid(m0)
      mr <- lmer(res ~ 1+(1|site),  data = data)
      summary(mr) # site explain 56% of var in residuals
         l=data.frame(summary(mr)$varcor)
         l = l[is.na(l$var2),]
         l$var1 = ifelse(is.na(l$var1),"",l$var1)
         l$pred = paste(l$grp,l$var1)
        ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
        ri$estimate_r = paste(ri$estimate_r,'%',sep='')
        ri    

# CHECK INDEPENDENCE OF RESIDUALS from Kubelka et al's original Table S6 models (Fig 3A)
  # TABLE S6a
    data <- read.csv( "DATApopulations.csv", sep=";")
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

    # Create a latitude variable

    data$NS <- data$Latitude > 0

    m0 <- lmekin( log(DPR) ~ (1|species) + log( N_nests) + mean_year * NS * abs( Latitude ), varlist = list( I, phyloMat, distanceMatrix), data = data )
      data$res = resid(m0)
        mr <- lmer(res ~ 1+(1|site),  data = data)
        summary(mr) # site explain 56% of var in residuals
           l=data.frame(summary(mr)$varcor)
           l = l[is.na(l$var2),]
           l$var1 = ifelse(is.na(l$var1),"",l$var1)
           l$pred = paste(l$grp,l$var1)
          ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
          ri$estimate_r = paste(ri$estimate_r,'%',sep='')
          ri  
    
  # TABLE S6b Historic
    data <- read.csv( "Before2000.csv", sep=";")
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

    # Create a latitude variable

    data$NS <- data$Latitude > 0


    hist(data$DPR_orig)
    # histogram of original daily predation rates, log transformation is really needed


    m0  <- lmekin( log(DPR) ~ (1|species) + log( N_nests) + NS + abs( Latitude ) + NS : abs( Latitude ), varlist = list( I, phyloMat, distanceMatrix), data = data )
        data$res = resid(m0)
        mr <- lmer(res ~ 1+(1|site),  data = data)
        summary(mr) # site explain 56% of var in residuals
           l=data.frame(summary(mr)$varcor)
           l = l[is.na(l$var2),]
           l$var1 = ifelse(is.na(l$var1),"",l$var1)
           l$pred = paste(l$grp,l$var1)
          ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
          ri$estimate_r = paste(ri$estimate_r,'%',sep='')
          ri  
    
    # TABLE S6c Recent
      data <- read.csv( "After2000.csv", sep=";")
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

      # Create a latitude variable

      data$NS <- data$Latitude > 0



      m0  <- lmekin( log(DPR) ~ (1|species) + log( N_nests) + NS + abs( Latitude ) + NS : abs( Latitude ), varlist = list( I, phyloMat, distanceMatrix), data = data )  
       data$res = resid(m0)
        mr <- lmer(res ~ 1+(1|site),  data = data)
        summary(mr) # site explain 56% of var in residuals
           l=data.frame(summary(mr)$varcor)
           l = l[is.na(l$var2),]
           l$var1 = ifelse(is.na(l$var1),"",l$var1)
           l$pred = paste(l$grp,l$var1)
          ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
          ri$estimate_r = paste(ri$estimate_r,'%',sep='')
          ri  

# try without too many years
  d = data
  dd =d[d$DPR_trans == 'NO' & d$Belt %in% c('Arctic','North temperate'),]
  dd =d[d$DPR_trans == 'NO' & d$Belt %in% c('Arctic','North temperate') & data$years_nr<10,]
  nrow(dd)
  m = lmer( log(DPR + 0.01) ~ log(N_nests) + mean_year  + Latitude + (1|site) , data = dd)
  m = lmer( log(DPR + 0.01) ~  mean_year  + (1|site) , data = dd)
  summary(m)
  summary(glht(m))

# Fig DPR original vs log
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
    png(glue(here::here(),"/Outputs/original_vs_ln-transformed_horizontal_base.png"), width=1.6*2,height=1.8,units="in",res=600) 
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
    mtext("ln-transformed scale", side = 3, line = -0.8, col='black', cex = 0.6, xpd = TRUE, at =mean(c(1944,2016)))
    
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

### NOT USED ###
# only for sites with multiple data
    tab <- table( data$site )
    idx <- which(tab > 1)
    inclPops <- names(tab)[idx]
    retain <- which( data$site %in% inclPops == TRUE)
    data_site_multiple <- data[retain,]

  mr2 <- lmer(res ~ 1+(1|site),  data = data_site_multiple)
  summary(mr2) # 72%
  l=data.frame(summary(mr2)$varcor)
     l = l[is.na(l$var2),]
     l$var1 = ifelse(is.na(l$var1),"",l$var1)
     l$pred = paste(l$grp,l$var1)
    ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
    ri$estimate_r = paste(ri$estimate_r,'%',sep='')
    ri

# only for sites with multiple data and no transformation      
  sub97 <- data[which(data$DPR_trans == "NO"),]
  tab <- table( sub97$site )
  idx <- which(tab > 1)
  inclPops <- names(tab)[idx]
  retain <- which( sub97$site %in% inclPops == TRUE)
  subsub97 <- sub97[retain,]

  mr2 <- lmer(res ~ 1+(1|site),  data = subsub97)
  summary(mr2) # 72%
  l=data.frame(summary(mr2)$varcor)
  	 l = l[is.na(l$var2),]
  	 l$var1 = ifelse(is.na(l$var1),"",l$var1)
  	 l$pred = paste(l$grp,l$var1)
  	ri=data.frame(type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
  	ri$estimate_r = paste(ri$estimate_r,'%',sep='')
  	ri


### RANDOM CRAP - to delete
      var_site = apply(bsim@ranef$site[,,1], 1, var)
      var_species = apply(bsim@ranef$species[,,1], 1, var)
      var_res = bsim@sigma^2
      quantile(var_site/(var_site+var_species+var_res), prob=c(0.025,0.5,0.975))

      0.37/var(log(dd$DPR))
      0.27^2/var(log(dd$DPR))
      2.770e-01/(2.770e-01+6.247e-09+7.219e-02)
      quantile(apply(bsim@ranef$site[,,1], 1, var), prob=c(0.025,0.5,0.975))
      quantile(apply(bsim@ranef$species[,,1], 1, var), prob=c(0.025,0.5,0.975))
      quantile(bsim@sigma^2, prob=c(0.025,0.5,0.975))

      

       l=data.frame(summary(m)$varcor)
       l = l[is.na(l$var2),]
       l$var1 = ifelse(is.na(l$var1),"",l$var1)
       l$pred = paste(l$grp,l$var1)

      rv = c(quantile(apply(bsim@ranef$site[,,1], 1, var), prob=c(0.5)),quantile(apply(bsim@ranef$species[,,1], 1, var), prob=c(0.5)), quantile(bsim@sigma^2, prob=c(0.5)))
      rlwr = c(quantile(apply(bsim@ranef$site[,,1], 1, var), prob=c(0.025)),quantile(apply(bsim@ranef$species[,,1], 1, var), prob=c(0.025)), quantile(bsim@sigma^2, prob=c(0.025)))
      rupr = c(quantile(apply(bsim@ranef$site[,,1], 1, var), prob=c(0.975)),quantile(apply(bsim@ranef$species[,,1], 1, var), prob=c(0.975)), quantile(bsim@sigma^2, prob=c(0.975)))
      data.frame(model=name,type='random',effect=l$pred,estimate=rv, lwr=rlwr, upr=rupr)
      
     
# get names of random part variables
l=data.frame(summary(m)$varcor)
l = l[is.na(l$var2),]
l$var1 = ifelse(is.na(l$var1),"",l$var1)
#l$pred = paste(l$grp,l$var1)

q050={}
q025={}
q975={}
pred={}

for (ran in names(bsim@ranef)) {
  ran_type = l$var1[l$grp == ran]
  for(i in ran_type){
    q050=c(q050,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.5)))
    q025=c(q025,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.025)))
    q975=c(q975,quantile(apply(bsim@ranef[[ran]][,,ran_type], 1, var), prob=c(0.975)))
    pred= c(pred,paste(ran, i))
  }
}
q050=c(q050,quantile(bsim@sigma^2, prob=c(0.5)))
q025=c(q025,quantile(bsim@sigma^2, prob=c(0.025)))
q975=c(q975,quantile(bsim@sigma^2, prob=c(0.975)))
pred= c(pred,'Residual')

name = 'bla'
data.frame(model=name,type='random',effect=pred,estimate=q050, lwr=q025, upr=q975)
xx = m_out(model = m)
m_out = function(name = "define", model = m, round_ = 3, nsim = 5000, aic = TRUE)
