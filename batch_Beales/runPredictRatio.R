setwd("D:/APAData/GitHub/River Fluxes/ratioestimator")

source('altmod.R')
source('collapse_stratbins.R')
source('mod.R')
source('nsamp.R')
source('predict_ratio.R')

#objects needed
siteQ<-read.csv("siteQ.csv")
siteConstit<-read.csv("siteConstit.csv")
siteQ$date<-as.Date(siteQ$date)
siteConstit$date<-as.Date(siteConstit$date)

minDaysPerYear<-345
constitName<-"NO3"
stationIdname<-"CODIGO_ESTACAO"
hi_flow_percentile <- 80 #Default threshold for designating high-flow observations */
ratio_strata_nsamp_threshold <- 10 #Default minimum number of observations required for inclusion of a stratum in the ratio estimate */
concTrans <- NA #constant transformation factor for converting to units of mg/L
qTrans<- 35.31466 #constant transformation factor for converting Q to units of ft3/s
loadUnits <- "kg"

predict_ratio(
  siteQ,siteConstit,minDaysPerYear,constitName,stationIdname,
  hi_flow_percentile,ratio_strata_nsamp_threshold,concTrans, qTrans)
