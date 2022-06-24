# LIBRARIES ----------------------------------------------------------
#Reading and Reporting
library(readxl) # read excel files
library(jtools) #for robust errors
library(sjPlot) #for tab_model

#Model comparison
library(performance) # model performance
library(lmtest) #Log likelihood ratio test
library (olsrr) #Collinearity Detection

# Install libraries for Spatial Regression 
library (sp) # dependecy for rgdal
library (rgdal) # read the shp file
library (terra) # dependecy for spdep
library(spdep) # for spatial weights generation
library (spatialreg)  # spatial regression
library (sphet) # spatial lag and error HET
library (spsur) #for sarar


# READ DATA -----------------------------------------------------------

# Working directory
setwd("C:/Users/Norbert/OneDrive/A08. Codruta si Stefana - Covid/03. Model")

#Sets the number of digits to 3, and scientific notation to 7
options(digits=3, scipen=3)

#Read data
model  <- read_excel("C:/Users/Norbert/OneDrive/A08. Codruta si Stefana - Covid/03. Model/data.xlsx", sheet = "data")
names(model)

# Spatial data preparation
uat = readOGR(dsn = ".", layer = "ro_atu") #citeste shp UAT-uri
names(uat) #afiseaza variabilele
uat$siruta=as.character(uat$siruta) # recodeaza SIRUTA in dimensiune

# Lists of spatial weights from the shape file DIFF
spatdata <- merge(uat, model, by = "siruta", all.x=F, all.y=F) #merge files
queen = poly2nb(spatdata) #neighbours list with single shared boundary point 
listw = nb2listw(queen) #converts neighbours list to listwise oject type

# CORRELATIONS --------------------------------------------------------------
# Calculate correlation matrixes
cormat <- cor(model[4:27])
t<-sjp.corr(model[4:27])

# Report correlation
t 
cormat


# EQUATIONS --------------------------------------------------------------
eq1 = VaccinationRateT1 ~ FamilyDoctors + Employees + ElderlyRetired + Migration + SocialistInv + RelativePoverty + ElderlyMinimum + Roma + LivingSpace2020 + PrimaryEd4 + MaxCovid + Neoprotestants + VotesYES2018 + Urban + Size100 + AirportsPublic
eq2 = VaccinationRateT1 ~ FamilyDoctors +                            + Migration + SocialistInv + RelativePoverty + ElderlyMinimum + Roma + LivingSpace2020 + PrimaryEd4 + MaxCovid + Neoprotestants + VotesYES2018 + Urban + Size100 + AirportsPublic
eq3 = VaccinationRateT2 ~ FamilyDoctors + Employees + ElderlyRetired + Migration + SocialistInv + RelativePoverty + ElderlyMinimum + Roma + LivingSpace2020 + PrimaryEd4 + MaxCovid + Neoprotestants + VotesYES2018 + Urban + Size100 + AirportsPublic
eq4 = VaccinationRateT2 ~ FamilyDoctors +                            + Migration + SocialistInv + RelativePoverty + ElderlyMinimum + Roma + LivingSpace2020 + PrimaryEd4 + MaxCovid + Neoprotestants + VotesYES2018 + Urban + Size100 + AirportsPublic

# OLS REGRESSION ---------------------------------------------------------

#Estimation

olsT1 <- lm(eq1, data=spatdata)
olsT2 <- lm(eq3, data=spatdata)

#Results

summ (olsT1, digits=3, robust=FALSE, vifs=TRUE)
summ (olsT2, digits=3, robust=FALSE, vifs=TRUE)
tab_model (olsT1, olsT2)

# Fit measures
compare_performance (olsT1, olsT2)
logLik(olsT1)
logLik(olsT2)
lrtest (olsT1, olsT2) 
ols_test_breusch_pagan(olsT1)
ols_test_breusch_pagan(olsT2)
ols_eigen_cindex(olsT1)
ols_eigen_cindex(olsT2)

#Moran test for olsT1 & olsT2
lm.morantest(olsT1,listw)
lm.morantest(olsT2,listw)

#Lagrange multiplier olsT1 & olsT2
lm.LMtests(olsT1,listw,test="all")
lm.LMtests(olsT2,listw,test="all")

#SARAR HET #------------------------------------------------------------
# Generalized Spatial two stage least squares estimation of a Cliff-Ord type model with Heteroskedastic Innovations
# Arraiz, I. and Drukker, M.D. and Kelejian, H.H. and Prucha, I.R. (2007)

# Model estimation
SararHetT1A=sphet::gstslshet(eq1, data=spatdata,listw)
SararHetT1B=sphet::gstslshet(eq2, data=spatdata,listw)
SararHetT2A=sphet::gstslshet(eq3, data=spatdata,listw)
SararHetT2B=sphet::gstslshet(eq4, data=spatdata,listw)


# Model report
summary (SararHetT1A, digits=3)
summary (SararHetT1B, digits=3)
summary (SararHetT2A, digits=3)
summary (SararHetT2B, digits=3)


# Model fit
SararHetT1A.R2<-1-sum(SararHetT1A$residuals^2)/(sum((spatdata$VaccinationRateT1-mean(spatdata$VaccinationRateT1))^2))
SararHetT1B.R2<-1-sum(SararHetT1B$residuals^2)/(sum((spatdata$VaccinationRateT1-mean(spatdata$VaccinationRateT1))^2))
SararHetT2A.R2<-1-sum(SararHetT2A$residuals^2)/(sum((spatdata$VaccinationRateT2-mean(spatdata$VaccinationRateT2))^2))
SararHetT2B.R2<-1-sum(SararHetT2B$residuals^2)/(sum((spatdata$VaccinationRateT2-mean(spatdata$VaccinationRateT2))^2))
SararHetT1A.R2
SararHetT1B.R2
SararHetT2A.R2
SararHetT2B.R2


# Salvare variabila prezisa
SararHet.pred <- data.frame(spatdata$siruta, SararHetT1A$yhat, SararHetT2A$yhat)
SararHet.name <- list("SIRUTA","YHat.T1", "YHat.T2")
colnames(SararHet.pred)[0] <- "ID"
write.csv (SararHet.pred, file="C:/Users/Norbert/OneDrive/A08. Codruta si Stefana - Covid/03. Model/SararHetT1.csv")
names(SararHet.pred)

#SUR-SARAR
#------------------------------------------------------------
# Seamingly Unrelated Regression model with spatial lags of the explained variables (with Spatial Autoregressive term) and spatial errors (and Spatial Autoregressive Disturbances)
# LeSage, J., and Pace, R. K. (2009)
# Mur, J., Lopez, F., and Herrera, M. (2010)

# T1 redus
SararHetT1.spsur <- spsur::spsurml(formula=eq2, listw=listw, type="sarar", data = spatdata)
summary(SararHetT1.spsur)

# library complet
SararHetT1.spsur.b <- spsur::spsurml(formula=eq1, listw=listw, type="sarar", data = spatdata)
summary(SararHetT1.spsur.b)

#save Yhat
model7A.spsur.pred <- data.frame (spat.2018$SIRUTA, model7A.spsur$fitted.values)
write.csv (model7A.spsur.pred, file="C:/Users/Norbert/OneDrive/A03. Cornel, Codruta - Depentent urbanization/02. Model/model7A.sarar.pred.csv")

