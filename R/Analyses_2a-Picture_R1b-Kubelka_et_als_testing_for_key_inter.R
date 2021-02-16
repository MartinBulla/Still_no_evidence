# ===========================================================================
# Supporting information for "Bulla, Valcu & Kempenaers. (2021) "Still no 
# evidence for disruption of global patterns of nest predation in shorebirds" 
# Contributor: Martin Bulla
# 📍 script runs relative to the project's root directory, uses Kubelka et 
# al's data and its manipulation, and generates Picture R1b within part
# (2) Testing the key interaction effect: does the temporal change in 
# predation vary across the globe? 
# ===========================================================================

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
  # [1] mgcv_1.8-31     nlme_3.1-148    ggplot2_3.3.2   phytools_0.7-47 maps_3.3.0      coxme_2.2-16   
  # [7] bdsmatrix_1.3-4 survival_3.1-12 ape_5.4-1      
  # 
  # loaded via a namespace (and not attached):
  #  [1] Rcpp_1.0.5              pillar_1.4.6            compiler_4.0.2          tools_4.0.2            
  #  [5] tibble_3.0.3            lifecycle_0.2.0         gtable_0.3.0            lattice_0.20-41        
  #  [9] pkgconfig_2.0.3         rlang_0.4.7             Matrix_1.2-18           fastmatch_1.1-0        
  # [13] igraph_1.2.5            parallel_4.0.2          expm_0.999-5            coda_0.19-3            
  # [17] withr_2.2.0             dplyr_1.0.1             generics_0.0.2          vctrs_0.3.2            
  # [21] gtools_3.8.2            tidyselect_1.1.0        combinat_0.0-8          grid_4.0.2             
  # [25] scatterplot3d_0.3-41    glue_1.4.2              R6_2.4.1                plotrix_3.7-8          
  # [29] animation_2.6           phangorn_2.5.5          purrr_0.3.4             magrittr_1.5           
  # [33] ellipsis_0.3.1          scales_1.1.1            MASS_7.3-51.6           splines_4.0.2          
  # [37] mnormt_2.0.2            colorspace_1.4-1        numDeriv_2016.8-1.1     quadprog_1.5-8         
  # [41] munsell_0.5.0           tmvnsim_1.0-2           crayon_1.3.4            clusterGeneration_1.3.4  # 