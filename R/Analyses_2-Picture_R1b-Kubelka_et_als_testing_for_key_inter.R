# TOOLS & SETTINGS from Kubelka et al.'s (2018) Models.R script from Dryad
  rm( list = ls() )

  library( ape )
  library( coxme )
  library( phytools )
  library(ggplot2)
  library(mgcv)

  # -----------------------------------------------
  # - Simple function for calculating distances ---
  # -----------------------------------------------
  geodetic <- function(l1, t1, l2, t2) {
    
    l1 <- l1 / 360.0 * 2.0 * pi
    l2 <- l2 / 360.0 * 2.0 * pi
    t1 <- t1 / 360.0 * 2.0 * pi
    t2 <- t2 / 360.0 * 2.0 * pi
    
    dist <- 6371.0 * acos( sin(t1) * sin(t2) + cos(t1) * cos(t2) * cos(l2 - l1) )
    
    return(dist)
    
  }

  # -----------------------------------------------
  # ------- Generates a distance matrix -------
  # -----------------------------------------------
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


  # Load up the data and trees
  trees <- read.tree( "Data/trees2.phy" )
  #data <- read.csv("Data/DATApopulations.csv", h = T, sep = ";")
  tree <- trees[[42]]  # just choose one to start with
  #plotTree(tree, node.numbers = TRUE, fsize = 0.5)

  # add in a tip:

  idx <- which(tree$tip == "Charadrius_alexandrinus")
  np <- nodepath(tree)[[idx]]
  to <- np[ length( np ) - 1] 
  ed <- which( tree$edge[,2] == idx)

  tree1 <- bind.tip( tree, "Charadrius_nivosus", edge.length = tree$edge.length[ed], where = to )
  #plotTree(tree1, node.numbers = TRUE, fsize = 0.5)


  addInTip <- function( phylo, where, newname) {
    idx <- which( phylo$tip == where )
    np <- nodepath( phylo)[[idx]]
    to <- np[ length(np) -1 ]
    ed <- which( tree$edge[,2] == idx )
    newphylo <- bind.tip( phylo, newname, edge.length = phylo$edge.length[ed], where = to )
    return( newphylo)
  }

  tree2 <- addInTip( tree1,"Gallinago_gallinago", "Gallinago_delicata" )
  #plotTree(tree2, node.numbers = TRUE, fsize = 0.5)

  ### data  
    data <- read.csv( "Data/DATApopulations.csv", sep=";")

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

# DPR model with interaction of interests and its unreported results (our Picture R1B) for DPR from Kubelka et al.'s (2018) Models.R script from Dryad
  # annotated as follows:
    ### -----------------------------------------------
    ### - Fig. 2A, 3A ; Table S6. Effect of latitude (A, B and C) and time (A) on nest predation --- daily nest predation rate - spatial and temporal patterns in one model
    ### -----------------------------------------------

    # all DPR values in the data are changed to original DPR + 0.01 to be able to use log(DPR) which is needed for data distribution normality.
  # models
    model_x <- lmekin( log(DPR) ~ (1|species) + log( N_nests) + mean_year * NS * abs( Latitude ), varlist = list( I, phyloMat, distanceMatrix), data = data )
    model_x # Picture R1B


    # next deleting nonsignificant interactions

    model_23 <- lmekin( log(DPR) ~ (1|species) + log( N_nests) + mean_year + NS + abs( Latitude ) + mean_year : NS + NS : abs( Latitude ) + mean_year : abs( Latitude ), varlist = list( I, phyloMat, distanceMatrix), data = data )
    model_23

# TPR model with interaction of interests and its unreported results from Kubelka et al.'s (2018) Models.R script from Dryad
  # annotated as follows: 
    ### -----------------------------------------------
    ### - Fig. 2B, 3B ; Table S6. Effect of latitude (A, B and C) and time (A) on nest predation --- total nest predation rate - spatial and temporal patterns in one model
    ### -----------------------------------------------

    # all TPR values in the data are original, no transformation

    model_y <- lmekin( TPR ~ (1|species) + log( N_nests) + mean_year * NS * abs( Latitude ), varlist = list( I, phyloMat, distanceMatrix), data = data )
    model_y


    # next deleting of nonsignificant interactions

    model_31 <- lmekin( TPR ~ (1|species) + log( N_nests) + mean_year + NS + abs( Latitude ) + mean_year : NS + NS : abs( Latitude ) + mean_year : abs( Latitude ), varlist = list( I, phyloMat, distanceMatrix), data = data )
    model_31

    par(mfrow = c(1, 2))
    plot(density(resid(model_31)))
    qqnorm(resid(model_31))
    qqline(resid(model_31))
    par(mfrow = c(1, 1))
    ## - check for data normality -> ok


    model_32 <- lmekin( TPR ~ (1|species) + log( N_nests) + mean_year + NS + abs( Latitude ) + NS : abs( Latitude ) + mean_year : abs( Latitude ), varlist = list( I, phyloMat, distanceMatrix), data = data )
    model_32

    par(mfrow = c(1, 2))
    plot(density(resid(model_32)))
    qqnorm(resid(model_32))
    qqline(resid(model_32))
    par(mfrow = c(1, 1))
    ## - check for data normality -> ok


    model_33 <- lmekin( TPR ~ (1|species) + log( N_nests) + mean_year + NS + abs( Latitude ) + mean_year : abs( Latitude ), varlist = list( I, phyloMat, distanceMatrix), data = data )
    model_33

    par(mfrow = c(1, 2))
    plot(density(resid(model_33)))
    qqnorm(resid(model_33))
    qqline(resid(model_33))
    par(mfrow = c(1, 1))
    ## - check for data normality -> ok


    # interaction are nonsignificant, now we will address main factors only

    model_34 <- lmekin( TPR ~ (1|species) + log( N_nests) + mean_year + NS + abs( Latitude ), varlist = list( I, phyloMat, distanceMatrix), data = data )
    model_34

    par(mfrow = c(1, 2))
    plot(density(resid(model_34)))
    qqnorm(resid(model_34))
    qqline(resid(model_34))
    par(mfrow = c(1, 1))
    ## - check for data normality -> ok