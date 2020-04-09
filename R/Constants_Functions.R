# CONSTANTS
  nsim = 5000
  recent = brewer.pal(n = 11, name = "Spectral")[2]#'#f6423f'#"#f04e47"#"#e95a54"#'#B2182B' #'#f74b45' ##'#f6423f' ##
  historic = "#2171B5"#"#2171B5"##"#6075b8"#"#5f72b1"#'#4393C3'# '#5b6ea3'
  rest = "#72c5c9"# "#909e85"
  rest_bg ="#f0f9f9"
  arct_bg = "#fceff5"
  col_p="gray53"  # color of point's outline
  col_pb="gray98"  # color of point's background	
  col_ = data.frame(
		year_ = c(1970,1980,1990,2000,2010),
		Belt = c('South temperate','South tropics','North tropics','North temperate','Arctic'),
		#line_col = c("lightseagreen","dimgray","brown","darkolivegreen","red"),
		line_col = c("#72c5c9","#9c9c9c","#bf786b","#729f44","#e7629c"),
		#line_col = c("lightseagreen","dimgray","brown","darkolivegreen","deeppink3"),
		poly_col = c("paleturquoise1","gainsboro","navajowhite2","darkolivegreen3","lightpink"),
		year_col = brewer.pal(n = 9, name = "Blues")[c(3,4,5,7,9)], #c('#e4d5de','#c5a6b9','#a67794','#7c516c','#4d3243'),#c('#D8C6D1','#CEB7C6','#BB9AAF','#A77E98','#946181'),#brewer.pal(n = 9, name = "Blues")[c(3,4,5,7,9)],#[c(5,6,7,8,9)],#
		#year_col = c('#5B6EA3','#866489','#B05B70','#DB5156','#F74B45'),#c('#D8C6D1','#CEB7C6','#BB9AAF','#A77E98','#946181'),
		#[c(5,6,7,8,9)],#
		stringsAsFactors=FALSE
		)
	
# FUNCTIONS from Kubelka et al.
  # Calculating distances 
    geodetic <- function(l1, t1, l2, t2) {
     l1 <- l1 / 360.0 * 2.0 * pi
     l2 <- l2 / 360.0 * 2.0 * pi
     t1 <- t1 / 360.0 * 2.0 * pi
     t2 <- t2 / 360.0 * 2.0 * pi
     dist <- 6371.0 * acos( sin(t1) * sin(t2) + cos(t1) * cos(t2) * cos(l2 - l1) )
     return(dist)
    }
  # Generating a distance matrix 
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
  # Adding species to the phylogenetic trees
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
    ed <- which( phylo$edge[,2] == idx ) # was phylo$tree and hence used only data from a signle tree
    newphylo <- bind.tip(tree = phylo, tip.label = newname, edge.length = phylo$edge.length[ed], where = to )
    return( newphylo)
  }
# model output function
  m_out = function(name = "define", model = m, round_ = 3, nsim = 5000, aic = TRUE){
	bsim <- sim(model, n.sim=nsim)  
	 v = apply(bsim@fixef, 2, quantile, prob=c(0.5))
	 ci = apply(bsim@fixef, 2, quantile, prob=c(0.025,0.975))	
	 oi=data.frame(model=name,type='fixed',effect=rownames(coef(summary(model))),estimate=v, lwr=ci[1,], upr=ci[2,])
			rownames(oi) = NULL
			oi$estimate_r=round(oi$estimate,round_)
			oi$lwr_r=round(oi$lwr,round_)
			oi$upr_r=round(oi$upr,round_)
	oii=oi[c('model','type',"effect", "estimate_r","lwr_r",'upr_r')]	
	
	 l=data.frame(summary(model)$varcor)
	 l = l[is.na(l$var2),]
	 l$var1 = ifelse(is.na(l$var1),"",l$var1)
	 l$pred = paste(l$grp,l$var1)
	ri=data.frame(model=name,type='random %',effect=l$pred, estimate_r=round(100*l$vcov/sum(l$vcov)), lwr_r=NA, upr_r=NA)
	ri$estimate_r = paste(ri$estimate_r,'%',sep='')
	x = rbind(oii,ri)
	if (aic == TRUE){	
		x$AIC = NA
		x$AIC[1]=AIC(model) # note that this one is incorrect because it is fitted with REML
		x$delta = x$prob = x$ER = NA
		}
    return(x)
  } 
# model assumption function
  m_ass = function(name = 'define', mo = m0, dat = d, fixed = NULL, categ = NULL, trans = NULL, spatial = TRUE, temporal = TRUE, PNG = TRUE){
   l=data.frame(summary(mo)$varcor)
   l = l[is.na(l$var2),]
   if(PNG == TRUE){
	png(paste(outdir,name, ".png", sep=""), width=6,height=9,units="in",res=600)
	 }else{dev.new(width=6,height=9)}
   
   n = nrow(l)-1+length(fixed)+length(categ) + 7
   par(mfrow=c(ceiling(n/3),3))
   
   scatter.smooth(fitted(mo),resid(mo),col='grey');abline(h=0, lty=2, col ='red')
   scatter.smooth(fitted(mo),sqrt(abs(resid(mo))), col='grey')
   qqnorm(resid(mo), main=list("Normal Q-Q Plot: residuals", cex=0.8),col='grey');qqline(resid(mo))
   #unique(l$grp[l$grp!="Residual"])
   for(i in unique(l$grp[l$grp!="Residual"])){
	#i = "mean_year"
	ll=ranef(mo)[names(ranef(mo))==i][[1]]
	if(ncol(ll)==1){
	 qqnorm(ll[,1], main = paste(i,names(ll)[1]),col='grey');qqline(ll[,1], col ='red')
	 }else{
	  qqnorm(ll[,1], main = paste(i,names(ll)[1]),col='grey');qqline(ll[,1], col ='red')
	  qqnorm(ll[,2], main = paste(i,names(ll)[2]),col='grey');qqline(ll[,2], col ='red')
	 }
	}
	
   # variables
   scatter={} 
   for (i in rownames(summary(mo)$coef)) {
		#i = "lat_abs"
      j=sub("\\).*", "", sub(".*\\(", "",i)) 
      scatter[length(scatter)+1]=j
    }
    x = data.frame(scatter=unique(scatter)[2:length(unique(scatter))],
					log_ = grepl("log",rownames(summary(mo)$coef)[2:length(unique(scatter))]), stringsAsFactors = FALSE)
    for (i in 1:length(fixed)){
	    jj =fixed[i]
		variable=dat[,names(dat)==jj]
		if(trans[i]=='log'){
        scatter.smooth(resid(mo)~log(variable),xlab=paste('log(',jj,')',sep=''), col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
		}else if(trans[i]=='abs'){
        scatter.smooth(resid(mo)~abs(variable),xlab=paste('abs(',jj,')',sep=''), col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
		}else{
        scatter.smooth(resid(mo)~variable,xlab=jj,col = 'grey');abline(h=0, lwd=1, lty = 2, col ='red')
      }
     }
	
	if(length(categ)>0){
	  for(i in categ){
		 variable=dat[,names(dat)==i]
		  boxplot(resid(mo)~variable, medcol='grey', whiskcol='grey', staplecol='grey', boxcol='grey', outcol='grey');abline(h=0, lty=3, lwd=1, col = 'red')
		 }
	}	  
		  
	if(temporal == TRUE){
		acf(resid(mo), type="p", main=list("Temporal autocorrelation:\npartial series residual",cex=0.8))
		}
	if(spatial == TRUE){	
	spdata=data.frame(resid=resid(mo), x=dat$Longitude, y=dat$Latitude)
		spdata$col=ifelse(spdata$resid<0,rgb(83,95,124,100, maxColorValue = 255),ifelse(spdata$resid>0,rgb(253,184,19,100, maxColorValue = 255), 'red'))
		#cex_=c(1,2,3,3.5,4)
		cex_=c(1,1.5,2,2.5,3)
		spdata$cex=as.character(cut(abs(spdata$resid), 5, labels=cex_))
	  plot(spdata$x, spdata$y,col=spdata$col, cex=as.numeric(spdata$cex), pch= 16, main=list('Spatial distribution of residuals', cex=0.8))
		legend("topleft", pch=16, legend=c('>0','<0'), ,col=c(rgb(83,95,124,100, maxColorValue = 255),rgb(253,184,19,100, maxColorValue = 255)), cex=0.8)
	  plot(spdata$x[spdata$resid<0], spdata$y[spdata$resid<0],col=spdata$col[spdata$resid<0], cex=as.numeric(spdata$cex[spdata$resid<0]), pch= 16, main=list('Spatial distribution of residuals (<0)', cex=0.8))
	  plot(spdata$x[spdata$resid>=0], spdata$y[spdata$resid>=0],col=spdata$col[spdata$resid>=0], cex=as.numeric(spdata$cex[spdata$resid>=0]), pch= 16, main=list('Spatial distribution of residuals (>=0)', cex=0.8))
		}
   
   mtext(paste(slot(mo,"call")[1],'(',slot(mo,"call")[2],sep=''), side = 3, line = -1, cex=0.7,outer = TRUE)
  if(PNG==TRUE){dev.off()}
  }
    