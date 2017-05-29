#### Setup for load modeling ####

library(loadflex)
library(EGRET)

meta <- metadata(
  constituent = 'NO3',
  consti.name = 'Nitrate',
  conc.units = 'mg L^-1', 
  flow = 'Q',
  flow.units = 'cms',
  load.units = 'kg', 
  load.rate.units = 'kg d^-1',
  dates = 'date',
  site.name = 'MOGU02900',
  site.id = 'MOGU02900')

# Create empty lists where we will store models and predictions
models <- list()


#### Load model input data ####

# Calibration/fitting data: date, concentration, and discharge. Used to fit a model
fitdat <- read.csv('three_ANA_sites/input/NO3/MOGU02900.csv')
fitdat$date <- as.Date(fitdat$date)
View(fitdat)
plotEGRET('plotConcTime', meta=meta, data=fitdat)
plotEGRET('plotConcQ', meta=meta, data=fitdat)
plotEGRET('plotFluxQ', meta=meta, data=fitdat)

# Estimation data: date and discharge. Used to generate daily predictions
estdat <- read.csv('three_ANA_sites/input/Q/MOGU02900.csv')
estdat$date <- as.Date(estdat$date)
View(estdat)
plot(log(Q) ~ date, estdat, type='l', main='Discharge vs Time', xlab='Date', ylab='ln(Q (cms))')
plotEGRET('multiPlotDataOverview', meta=meta, data=fitdat, newdata=estdat)


#### Types of load models ####

# Interpolation
models$INT <- loadInterp(
  data=fitdat, metadata=meta,
  interp.format='conc', interp.function=rectangularInterpolation)
plotEGRET('plotConcTimeDaily', load.model=models$INT, newdata=estdat); title('INT', line=-2)

# Regression (LOADEST 5 parameter model, L5)
models$REG <- loadReg2(
  loadReg(NO3 ~ model(7), data=fitdat, 
          flow='Q', dates='date', flow.units='cms', conc.units='mg/L', load.units='kg',
          station=getInfo(meta, 'site.name')),
  site.id=getInfo(meta, 'site.id'), consti.name='Nitrate', pred.format='conc')
plotEGRET('plotConcTimeDaily', load.model=models$REG, newdata=estdat); title('REG', line=-2)

# Composite (regression with correction)
models$CMP <- loadComp(
  reg.model=models$REG, interp.data=fitdat,
  store=c('data','fitting.function'), interp.format='conc',
  interp.function=rectangularInterpolation)
plotEGRET('plotConcTimeDaily', load.model=models$CMP, newdata=estdat); title('CMP', line=-2)

# Weighted regression on time, discharge, and season (WRTDS)
input_WRTDS <- convertToEGRET(data=fitdat, newdata=estdat, meta=meta)
models$WRTDS <- modelEstimation(input_WRTDS, minNumObs=30)
plotConcTimeDaily(models$WRTDS); title('WRTDS', line=-2)
EGRET::plotConcHist(models$WRTDS, plotFlowNorm=FALSE)
EGRET::plotConcHist(models$WRTDS, plotFlowNorm=TRUE)
EGRET::plotFluxHist(models$WRTDS, plotFlowNorm=FALSE)
EGRET::plotFluxHist(models$WRTDS, plotFlowNorm=TRUE)
