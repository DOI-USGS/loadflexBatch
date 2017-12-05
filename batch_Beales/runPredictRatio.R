bealesFiles <- sapply(file.path('batch_Beales', c('altmod.R','collapse_stratbins.R','mod.R','nsamp.R','predict_ratio.R')), source)

#objects needed
siteQ<-read.csv("D:/APAData/GitHub/River Fluxes/ratioestimator/siteQ.csv")
siteConstit<-read.csv("D:/APAData/GitHub/River Fluxes/ratioestimator/siteConstit.csv")
siteQ$date<-as.Date(siteQ$date)
siteConstit$date<-as.Date(siteConstit$date)

minDaysPerYear<-345
dateName<-"date"
qName<-"Q"
constitName<-"NO3"
stationIdname<-"CODIGO_ESTACAO"
hi_flow_percentile <- 80 #Default threshold for designating high-flow observations */
ratio_strata_nsamp_threshold <- 10 #Default minimum number of observations required for inclusion of a stratum in the ratio estimate */
concTrans <- NA #constant transformation factor for converting to units of mg/L
qTrans<- 35.31466 #constant transformation factor for converting Q to units of ft3/s
loadUnits <- "kg"

predict_ratio(
  siteQ,siteConstit,minDaysPerYear,dateName,qName,constitName,stationIdname,
  hi_flow_percentile,ratio_strata_nsamp_threshold,concTrans, qTrans)
