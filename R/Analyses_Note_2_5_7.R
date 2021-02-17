# ===========================================================================
# Supporting information for "Bulla, Valcu & Kempenaers. (2021) "Still no 
# evidence for disruption of global patterns of nest predation in shorebirds" 
# Contributor: Martin Bulla
# 📍 script runs relative to the project's root directory and uses Kubelka et 
# al's data to generate Note 2, 5, 7 and 9 results 
# ===========================================================================

# TOOLS and DATA manipulation (mostly) from Kubelka et al
  rm( list = ls() )
  source('R/Constants_Functions.R')
  # load additional package 
    library( ape )
    library( coxme )
    library( phytools )
    library(mgcv)
    require(effects)
    require(plyr)
    require(glue)
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
    #plotTree(tree1, node.numbers = TRUE, fsize = 0.5)
    tree2 <- addInTip( tree1,"Gallinago_gallinago", "Gallinago_delicata" )
  # load and prepare data
    data <- read.csv( "Data/DATApopulations.csv", h = T, sep=";")	
	  datapred <- read.csv("Kubelka et al. 2019 Response_Data, R codes and supporting information/DATApopulations.csv", sep =";", h = T)

    # prepare subsets
      datapred$pop_ID <-paste(datapred$Longitude, datapred$Latitude)
      datapred$site <-paste(datapred$Longitude, datapred$Latitude)
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

# Note 2 - scrutinizing Kubelka et al's Response models about temporal changes in predation
  # M2 - model used in Kubelka et al. (2019a) as a replication of our model (lmer( log(DPR + 0.01) ~ log(N_nests) + Latitude + mean_year + (1|pop_ID))
    m2 <- lmer( log(DPR_orig + 0.01) ~ mean_year + (1|site) , data = subsub97) # uses only data with multiple data points per site
    summary(m2)
    summary(glht(m2))
  
  # Table A - M2 fitted to all directly calculated estimates (as we did in our Comment)
    m2c_all <- lmer( log(DPR_orig + 0.01) ~  mean_year  + (1|site) , data = sub97) 
    summary(m2c_all)
    summary(glht(m2c_all))
    A = m_out(name = "A", model = m2c_all, round_ = 3, nsim = 5000, aic = FALSE )[1:6]    
    A$p_value = NA
    A$p_value[1:2] = round(summary(glht(m2c_all))$test$pvalues,2)
    A$N = NA
    A$N[1] = nrow(sub97)
    
  # Table B - M2 controlled for number of nests to generate an estimate and for latitude of the population
    m2c <- lmer( log(DPR_orig + 0.01) ~ log(N_nests)  + mean_year + Latitude + (1|site) , data = subsub97) 
    summary(m2c)
    summary(glht(m2c))
    B = m_out(name = "B", model = m2c, round_ = 3, nsim = 5000, aic = FALSE )[1:6]    
    B$p_value = NA
    B$p_value[1:4] = round(summary(glht(m2c))$test$pvalues,2)
    B$N = NA
    B$N[1] = nrow(subsub97)

  # Table C - M2 controlled for number of nests to generate an estimate and including interaction with latitude as in our original Figure 1F (referred to by Kubelka et al. 2019) 
    m2c_asF3E <- lmer( log(DPR_orig + 0.01) ~ log(N_nests) + scale(mean_year)*poly(Latitude,3) + (1|pop_ID) , data = subsub97) 
    summary(m2c_asF3E)
    summary(glht(m2c_asF3E))  

    C = m_out(name = "C", model = m2c_asF3E, round_ = 3, nsim = 5000, aic = FALSE )[1:6]    
    C$p_value = NA
    C$p_value[1:9] = round(summary(glht(m2c_asF3E))$test$pvalues,2)
    C$N = NA
    C$N[1] = nrow(subsub97)
  
  # Table D - M2 controlled for number of nests to generate an estimate and including interaction with latitude as in our original Figure 1F (referred to by Kubelka et al. 2019), but using lat as continuous
    m2c_asF3Elat <- lmer( log(DPR_orig + 0.01) ~ log(N_nests) + scale(mean_year)*scale(Latitude) + (1|pop_ID) , data = subsub97) 
    summary(m2c_asF3Elat)
    summary(glht(m2c_asF3Elat))

    D = m_out(name = "D", model = m2c_asF3Elat, round_ = 3, nsim = 5000, aic = FALSE )[1:6]    
    D$p_value = NA
    D$p_value[1:5] = round(summary(glht(m2c_asF3Elat))$test$pvalues,2)
    D$N = NA
    D$N[1] = nrow(subsub97)

   # export
    write_xlsx(rbind(A,B,C,D), 'Outputs/Table_Note2.xlsx')

# Note 5 - Figure: DPR original vs log scale
  d = data
  d$Belt = as.factor(d$Belt) # define site
  d$site = paste(d$Latitude,d$Longitude) # define site
  d$lat_abs = abs(d$Latitude) # abs latitude
  d$ln_N_nests = log(d$N_nests)
  d$hemisphere =as.factor(ifelse(d$Latitude > 0, "Northern", "Southern"))
  d$genus = gsub("\\_.*","",d$species)
  d$DPR_trans[which(d$DPR_trans == 'YES' & d$source_id=="209")] = "NO"# source_id 209 Schekkerman et al. 1998 (Calidris 

  # predictions for transformed
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
        newD$pred <- exp(X%*%v)-0.01 # #newD$fit_b <- plogis(X%*%v) # in case on binomial scaleback
        predmatrix <- matrix(nrow=nrow(newD), ncol=nsim)
        for(i in 1:nsim) predmatrix[,i] <- exp(X%*%bsim@fixef[i,])-0.01
            newD$lwr <- apply(predmatrix, 1, quantile, prob=0.025)
            newD$upr <- apply(predmatrix, 1, quantile, prob=0.975)
        newD$line_col = col_$line_col[match(newD$Belt, col_$Belt)]
    dprC=newD   
  
  # predictions for original
    dd =d[d$DPR_trans == 'NO' & d$Belt %in% c('Arctic','North temperate'),]
    nrow(dd)
    #summary(dd$Belt) 
    m = lmer(DPR_orig ~ ln_N_nests + poly(mean_year,2)*Belt+(1|site)+(1|species),  data =dd )
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
    png(glue(here::here(),"/Outputs/Figure_Note5_original_vs_log-transformed.png"), width=1.6*2,height=1.8,units="in",res=600) 
     }else{dev.new(width=1.6*2,height=1.8)}
    
    par(oma = c(1.0, 1.5, 0.55, 0.25),mgp=c(1.2,0.15,0), 
        ps=12,font.main = 1, las=1, lwd = 0.5, tcl=-0.05,
        cex=1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5,
        col.lab="black", col.main="162", fg="grey45",
        bty="n",xpd=TRUE) #
  
    layout(mat = matrix(c(1,2),ncol = 2))
     
    pp = dprCorig

    par(mar=c(0.1,0.5,0,0.1),ps=12, cex=1, font.main = 1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.05,bty="n",xpd=FALSE)
   
    plot(pp$pred~pp$mean_year, pch=19,xlim=c(1944,2016), ylim=c(0,0.06),ylab=NA,xlab=NA, xaxt='n',yaxt='n', type='n')#xaxs="i",yaxs="i")
    
    axis(1, at=c(1944,1960,1980,2000,2016),labels=c(1944,1960,1980,2000,2016),cex.axis=0.5,mgp=c(0,-0.35,0), lwd = 0.5)
    axis(2, at=seq(0,0.06, by = 0.01), labels=c('0','0.01','0.02','0.03','0.04','0.05','0.06'), lwd = 0.5)
    mtext('Daily nest predation',side=2,line=1.1, cex=0.6, las=3, col = 'black')
    #text(x = 2016+5, y =0.07*1.01, labels= expression(bold("Based on model with predation on")), col='black', cex = 0.6, xpd = TRUE)
    mtext('Year',side=1,line=0.2, cex=0.6, las=1, col = 'black', at = 2016+10 )
    mtext(expression(bold("Based on model with predation rate on:")), side = 3, line = -0.2, col='black', cex = 0.6, xpd = TRUE, at = 2016+10)
    mtext("original-scale", side = 3, line = -0.8, col='black', cex = 0.6, xpd = TRUE, at =mean(c(1944,2016)))
    #text(x = 2016-30, y =0.07*0.98, labels= expression(bold("original scale")), col='black', cex = 0.6, xpd = TRUE)
    for(i in unique(pp$Belt)){
    #i ="Arctic"
    print(i)
    ppi = pp[pp$Belt==i,]
    polygon(c(ppi$mean_year, rev(ppi$mean_year)), c((ppi$lwr)-0.01, 
      rev((ppi$upr)-0.01)), border=NA, col = adjustcolor(ppi$line_col,alpha.f = 0.1))#adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
    lines(ppi$mean_year, (ppi$pred)-0.01, col=ppi$line_col,lwd=1)
            } 

    text(x =1990, y = 0.035, labels = "Arctic", cex = 0.5, col = col_$line_col[col_$Belt == "Arctic"])        
    text(x =1990, y = 0.0225, labels = "North temperate", cex = 0.5, col = col_$line_col[col_$Belt == "North temperate"])        
    
    pp = dprC
    #par(mar=c(0.1,0.5,0,0.1),ps=12, cex=1, font.main = 1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.05,bty="n",xpd=FALSE)
    par(mar=c(0.1,0.4,0,0.2),ps=12, cex=1, font.main = 1, cex.lab=0.6,cex.main=0.7, cex.axis=0.5, tcl=-0.05,bty="n",xpd=FALSE)
   
    plot(pp$pred~pp$mean_year, pch=19,xlim=c(1944,2016), ylim=c(0,0.06),ylab=NA,xlab=NA, xaxt='n',yaxt='n', type='n')#xaxs="i",yaxs="i")
    
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
      polygon(c(ppi$mean_year, rev(ppi$mean_year)), c(ppi$lwr, 
        rev(ppi$upr)), border=NA, col = adjustcolor(ppi$line_col,alpha.f = 0.1))#adjustcolor(col_t ,alpha.f = 0.2)) #0,0,0 black 0.5 is transparents RED
      lines(ppi$mean_year, ppi$pred, col=ppi$line_col,lwd=1)
              } 
    if(PNG == TRUE){dev.off()}
  
    #quartz.save(glue(here::here(),"/Outputs/original_vs_ln-transformed_horizontal_quartz.png"), type = "png") 
  
# Note 7 - scrutinizing Kubelka et al's response models about general increase in predation rates over time
  # for Kubelka et al's Table S3 results, but controlled for site see script Analyses_3a  

  # M3 not controlled for multiple sampling per site
      M3 <- lm( log( DPR) ~ log(N_nests) + mean_year +  Latitude , data = sub97) 
      summary(M3)
  # M4 adds the controls for multiple sampling per site
      M4 <- lmer( log(DPR) ~ log(N_nests) + mean_year  +  Latitude + (1|site) , data = sub97) 
      summary(M4)

      A = m_out(name = "A", model = M4, round_ = 3, nsim = 5000, aic = FALSE )[1:6]    
      A$p_value = NA
      A$p_value[1:4] = round(summary(glht(M4))$test$pvalues,2)
      A$N = NA
      A$N[1] = nrow(sub97)
      
      write_xlsx(A, 'Outputs/Table_Note7.xlsx')
  
  # the key is the control for site, as year remains non-significant when latitude is removed and even when log(N_nests) is removed
        m <- lmer( log(DPR) ~ log(N_nests) + mean_year  + (1|site) , data = sub97) # simplified
        m<- lmer( log(DPR) ~ mean_year  + (1|site) , data = sub97)  #M4 not controlled for number of nests for a given estimate and latitude 
        summary(glht(m))
        summary(M4)


        tab <- table( sub97$pop_ID )
        idx <- which(tab > 1)
        inclPops <- names(tab)[idx]


        retain <- which( sub97$pop_ID %in% inclPops == TRUE)
        subsub97 <- sub97[retain,]
        model.lmer <- lmer( log(DPR) ~ mean_year + (1|pop_ID) , data = subsub97) 
        summary(glht(model.lmer))

        m <- lmer( log(DPR) ~ log(N_nests) + Latitude + mean_year + (1|pop_ID) , data = subsub97) 
        summary(glht(m))
        m <- lmer( log(DPR) ~ log(N_nests) + mean_year  + (1|pop_ID) , data = subsub97) 
        m <- lmer( log(DPR) ~ Latitude + mean_year  + (1|pop_ID) , data = subsub97) 

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
  #  [1] lmerTest_3.1-3     glue_1.4.2         plyr_1.8.6         effects_4.1-4      carData_3.0-4     
  #  [6] mgcv_1.8-31        nlme_3.1-148       phytools_0.7-47    maps_3.3.0         coxme_2.2-16      
  # [11] bdsmatrix_1.3-4    ape_5.4-1          writexl_1.3.1      viridis_0.5.1      viridisLite_0.3.0 
  # [16] RColorBrewer_1.1-2 performance_0.4.8  multcomp_1.4-13    TH.data_1.0-10     survival_3.1-12   
  # [21] mvtnorm_1.1-1      ggplot2_3.3.2      data.table_1.13.0  arm_1.11-2         lme4_1.1-23       
  # [26] Matrix_1.2-18      MASS_7.3-51.6     
  # 
  # loaded via a namespace (and not attached):
  #  [1] insight_0.9.0           numDeriv_2016.8-1.1     tools_4.0.2             backports_1.1.8        
  #  [5] R6_2.4.1                rpart_4.1-15            DBI_1.1.0               Hmisc_4.4-0            
  #  [9] colorspace_1.4-1        nnet_7.3-14             withr_2.2.0             tidyselect_1.1.0       
  # [13] gridExtra_2.3           mnormt_2.0.2            phangorn_2.5.5          compiler_4.0.2         
  # [17] htmlTable_2.0.1         animation_2.6           expm_0.999-5            sandwich_2.5-1         
  # [21] bayestestR_0.7.2        scales_1.1.1            checkmate_2.0.0         quadprog_1.5-8         
  # [25] stringr_1.4.0           digest_0.6.25           foreign_0.8-80          minqa_1.2.4            
  # [29] base64enc_0.1-3         jpeg_0.1-8.1            pkgconfig_2.0.3         htmltools_0.5.0        
  # [33] plotrix_3.7-8           htmlwidgets_1.5.1       rlang_0.4.7             rstudioapi_0.11        
  # [37] generics_0.0.2          zoo_1.8-8               combinat_0.0-8          gtools_3.8.2           
  # [41] acepack_1.4.1           dplyr_1.0.1             magrittr_1.5            Formula_1.2-3          
  # [45] Rcpp_1.0.5              munsell_0.5.0           abind_1.4-5             lifecycle_0.2.0        
  # [49] scatterplot3d_0.3-41    stringi_1.5.3           clusterGeneration_1.3.4 grid_4.0.2             
  # [53] parallel_4.0.2          crayon_1.3.4            lattice_0.20-41         splines_4.0.2          
  # [57] tmvnsim_1.0-2           knitr_1.29              pillar_1.4.6            igraph_1.2.5           
  # [61] boot_1.3-25             codetools_0.2-16        fastmatch_1.1-0         mitools_2.4            
  # [65] latticeExtra_0.6-29     png_0.1-7               vctrs_0.3.2             nloptr_1.2.2.2         
  # [69] gtable_0.3.0            purrr_0.3.4             xfun_0.16               survey_4.0             
  # [73] coda_0.19-3             tibble_3.0.3            cluster_2.1.0           statmod_1.4.34         
  # [77] ellipsis_0.3.1 