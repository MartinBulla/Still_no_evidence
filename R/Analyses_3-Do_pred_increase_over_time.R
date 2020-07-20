# load packages, constants and data
    rm( list = ls() )   
    source('R/Constants_Functions.R')
    d = read.csv("Data/DATApopulations.csv", h = T, sep=";",stringsAsFactors = FALSE)
    d$site = paste(d$Latitude,d$Longitude) # define site
    dn = d[d$DPR_tran == 'NO',] # only high quality data
  
# Kubelka et al.'s (2018) Table S3a model	
     mk = lm(log(DPR) ~ log( N_nests) + mean_year,  data = d)
     summary(mk)
# Kubelka et al.'s (2018) Table S3a model CONTROLLED FOR STUDY SITE - effect reduced by 30%
     mkc = lmer(log(DPR) ~ log( N_nests) + mean_year+(1|site),  data = d) 
     summary(mkc)

# Kubelka et al.'s (2018) Table S3d model on non-transformed data 
     mkn = lm(log(DPR) ~ log( N_nests) +mean_year,  data = dn)  
     summary(mkn)
# Kubelka et al.'s (2018) Table S3d model on non-transformed data CONTROLLED FOR STUDY SITE - effect reduced by 30%
     mknc = lmer(log(DPR) ~ log( N_nests) +mean_year   +(1|site),  data = dn)  
     summary(mknc)
     m_out(name = "kmnc", model = mknc, round_ = 3, nsim = 5000, aic = FALSE)    