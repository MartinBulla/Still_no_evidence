# Global predation rates increase over time?


rm(list = ls())

library(lme4)
library(lmerTest)
datapred <- read.csv("Kubelka et al. 2019 Response_Data, R codes and supporting information/DATApopulations.csv", sep =";", h = T)


datapred$pop_ID <-paste(datapred$Longitude, datapred$Latitude)
datapred$hemis <-  datapred$Latitude > 0 
datapred$arctic <- datapred$Belt == "Arctic"
sub97 <- datapred[which(datapred$DPR_trans == "NO"),]

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
summary(model.lmer)

m <- lmer( log(DPR + 0.01) ~ log(N_nests) + Latitude + mean_year + (1|pop_ID) , data = subsub97) 
summary(glht(m))
m <- lmer( log(DPR + 0.01) ~ log(N_nests) + mean_year  + (1|pop_ID) , data = subsub97) 
m <- lmer( log(DPR + 0.01) ~ Latitude + mean_year  + (1|pop_ID) , data = subsub97) 




