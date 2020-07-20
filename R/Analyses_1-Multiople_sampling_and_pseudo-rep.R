# (1) Multiple sampling and pseudo-replication

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

# Figure R1 - distribution of sites over globe and sample sizes 
  xx = ddply(data, .(Latitude, Longitude, period_orig), summarise, N_estimates = length(Latitude))
  ggplot(xx, aes(y = Latitude, x = Longitude, col = period_orig)) + 
      geom_point(aes(size = N_estimates), alpha = 0.5) + 
      scale_size(name = "# of estimates per site", breaks = c(1, 5, 10), limits = c(1, 15))+ 
      scale_color_manual(values = c('red','blue'), name = "Mean year of the study", labels =c("\u2265 2000", "< 2000")) +
      guides(size = guide_legend(order = 1),col = guide_legend(order = 2))
  ggsave(file=glue(here::here(),"/Outputs/N-estimates_per_site_geography_period.png"), dpi = 300, width = 12, height = 8, units = "cm")

  # other versions
    xx = ddply(data, .(Latitude, Longitude), summarise, N_estimates = length(Latitude))
    xy = ddply(data, .(Latitude, Longitude, mean_year, period_orig, Belt), summarise, N_estimates = length(Latitude), med_TPR = median(TPR), avg_TPR = mean(TPR), med_DPR = median(DPR), avg_DPR = mean(DPR))
    
    ggplot(data, aes(y = Latitude, x = Longitude)) + geom_jitter(aes(size = N_nests), alpha = 0.5)
    quartz.save(glue(here::here(),"/Outputs/N-nest_per_estimate_geography.png"), type = "png")
    ggsave(file = glue(here::here(),"/Outputs/N-nest_per_estimate_geography.png", dpi = 300, width = 10, height = 8, units = "cm") 
    
    ggplot(xx, aes(y = Latitude, x = Longitude)) + geom_point(aes(size = N_estimates), alpha = 0.5) + scale_size(limits = c(1, 15)) + scale_size(name = "# of estimates")
    quartz.save(glue(here::here(),"/Outputs/N-estimates_per_site_geography.png"), type = "png")
    ggsave(file=glue(here::here(),"/Outputs/N-estimates_per_site_geography.png"), dpi = 300, width = 10, height = 8, units = "cm")
    

    ggplot(data, aes(x = abs(Latitude), y = TPR, col = mean_year)) + geom_jitter(aes(size = N_nests), alpha = 0.5)
    ggplot(data, aes(x = abs(Latitude), y = DPR, col = mean_year)) + geom_jitter(aes(size = N_nests), alpha = 0.5)
   
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

# Table R1 and further data summaries
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