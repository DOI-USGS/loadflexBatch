# This script runs loadflex in "batch mode": It fits several load estimation 
# models for many constituents at many sites, generates and saves the 
# predictions, and produces summaries over all models, constituents, and sites. 
# See https://github.com/USGS-R/loadflexBatch/blob/master/blog.md for an
# overview and instructions on how to use this file.

#------------------User Inputs--------------------#

inputFolder <- "three_ANA_sites/input" #folder containing all input subfolders

#input constituents
#script will look for folders with this name inside inputFolder, 
#and use in site metadata
constituents <- c("NO3", "PT")
loadUnits <- "kg"
loadRateUnits <- "kg/d"

dischargeFolder <- "Q" #subfolder of inputFolder containing discharge measurements for predictions
siteInfo <- "siteInfo.csv" #also inside inputFolder, data frame of site info

outputFolder <- "three_ANA_sites/output"  #output files and subfolders created here

#-------------------------Load packages, check files, set up directories-----------------------# 

library(dplyr)
library(loadflex)
library(tools)
library(rloadest)
source('batchHelperFunctions.R') #functions stored here

#need at least rloadest 0.4.4 for formula fix
if(compareVersion(as.character(packageVersion("rloadest")),"0.4.4") == -1) {
  stop("rloadest version 0.4.4 or greater is required")
}

fileDF <- makeFileDF(inputFolder, constits = constituents, discharge.folder = dischargeFolder)
allSiteInfo <- read.csv(file.path(inputFolder, siteInfo), stringsAsFactors = FALSE)

#setup output directories
nConstits <- length(constituents)
outConstit <- file.path(rep(outputFolder, nConstits), constituents)
sapply(outConstit, dir.create, recursive = TRUE, showWarnings = FALSE)
outTemporal <- file.path(rep(outConstit,4), c(rep("inputs", nConstits), 
                                              rep("annual", nConstits), 
                                              rep("multiYear", nConstits), 
                                              rep("modelMetrics", nConstits)))
sapply(outTemporal, dir.create, showWarnings = FALSE)

lastConstit <- NULL
graphics.off() #don't want open PDF connections

#-----------------loadflex--------------#

allModels <- list() # this might get big. i'd prefer splitting & saving

#loop over unique sites
for(i in 1:nrow(fileDF)) {
  message(paste('processing constituent file', fileDF$constitFile[i], '\n'))
  
  #TODO: use siteInfo to get column names, fail if they aren't the same as directories
  
  #read in appropriate files
  siteQ <- read.csv(fileDF$qFile[i], stringsAsFactors = FALSE)
  siteConstit <- read.csv(fileDF$constitFile[i], stringsAsFactors = FALSE)
  
  #convert date text to Dates 
  siteConstit$date <- as.Date(siteConstit$date)
  siteQ$date <- as.Date(siteQ$date)
  
  #pull out appropriate rows of allSiteInfo for Q and constit
  #need to extract constit and sites from file paths, so we know 
  #what row of site info to look at
  constitSite <- basename(file_path_sans_ext(fileDF$constitFile[i])) 
  constitName <- basename(dirname(fileDF$constitFile[i]))
 
  #if switching to a new consituent, open a new pdf
  #important that fileDF is sorted by consituent!
  #this should be the case the way makeFileDF looks at folders
  if(is.null(lastConstit)) {
    pdf(height = 11, width = 8.5, 
        file = file.path(outputFolder, constitName, sprintf("%s_plots.pdf", constitName)))
    lastConstit <- constitName
  }
  if(constitName != lastConstit) {
    lastConstit <- constitName
    dev.off()
    pdf(height = 11, width = 8.5, 
        file = file.path(outputFolder, constitName, sprintf("%s_plots.pdf", constitName)))
  }
  
  constitSiteInfo <- filter(allSiteInfo, matching.site == constitSite, constituent == constitName)
  qSiteInfo <- filter(allSiteInfo, matching.site == constitSite, constituent == 'Q')
  
  #deal with different discharge/consituent drainage areas
  if(qSiteInfo$basin.area != constitSiteInfo$basin.area) {
    ratio <- constitSiteInfo$basin.area/qSiteInfo$basin.area
    #modify discharge for both DFs
    siteConstit$Q <- siteConstit$Q*ratio
    siteQ$Q <- siteQ$Q*ratio
    message("Scaling discharge by basin area")
  }
  
  #create metadata
  #not sure units etc are following the correct format
  #currently depending on consistent column positions to get names
  constitColName <- names(siteConstit)[3]
  qwconstitColName <- paste0(constitColName, '_qw')
  qColName <- names(siteConstit)[2]
  dateColName <- names(siteConstit)[1]
  
  # format censored data for rloadest. For ANA, Status 0 means null or blank. 
  # Status 1 means a valid value and 2 means that the respective value is a 
  # detection limit."
  censor.statuses = c('0'='', '1'='', '2'='<')
  siteConstit[[qwconstitColName]] <- 
    smwrQW::as.lcens(
      values=ifelse(siteConstit[['status']] == 0, NA, siteConstit[[constitColName]])/2, 
      detlim=ifelse(siteConstit[['status']] < 2, 0, siteConstit[[constitColName]]), 
      censor.codes=censor.statuses[as.character(siteConstit[['status']])])
  siteConstit[[constitColName]] <- ifelse(
    siteConstit[['status']] == 0, NA, 
    ifelse(siteConstit[['status']] == 2, siteConstit[[constitColName]] / 2, # it's still bad, but use half MDL because better than full MDL
           siteConstit[[constitColName]]))
  # remove NA values, some of which may have been added by checking status
  siteConstit <- siteConstit[!is.na(siteConstit[[constitColName]]), ]
  
  # create a formal metadata object. site.id and flow.site.id must both equal
  # constitSite for our input file scheme to work
  siteMeta <- metadata(
    constituent = constitColName, consti.name = constitColName, conc.units = constitSiteInfo$units, 
    flow = qColName, flow.units = qSiteInfo$units, 
    load.units = loadUnits, load.rate.units = loadRateUnits, dates = dateColName,
    site.name = constitSiteInfo$site.name, site.id = constitSiteInfo$site.id, lat = constitSiteInfo$lat, lon = constitSiteInfo$lon, basin.area = constitSiteInfo$basin.area,
    flow.site.name = qSiteInfo$site.name, flow.site.id = qSiteInfo$site.id, flow.lat = qSiteInfo$lat, flow.lon = qSiteInfo$lon, flow.basin.area = qSiteInfo$basin.area
  )
  
  # compute and save info on the site, constituent, and input datasets (we'll 
  # recombine in the next loop). compute num.censored specially here because 
  # we're using rloadest format for censored data (smwrQW format) and will 
  # eventually have something simpler in place for loadflex, at which point
  # we'll add that to summarizeInputs.
  inputMetrics <- summarizeInputs(siteMeta, fitdat=siteConstit, estdat=siteQ)
  inputMetrics$fitdat.num.censored <- length(which(!is.na(siteConstit[[qwconstitColName]]@.Data[,'detlim'])))
  inputMetrics$estdat.num.censored <- NULL # assuming there isn't and shouldn't be censoring in Q. is that right?
  write.csv(inputMetrics, file.path(outputFolder, constitName, "inputs", paste0(constitSite, '.csv')), row.names=FALSE)
  
  #fit models
  set.seed(9451)
  #TODO: decide on standard column names?  user input timestep above?
  siteConstitRloadest <- setNames(
    siteConstit[c(dateColName, qwconstitColName, qColName)],
    c(dateColName, constitColName, qColName))
  loadRegFormula <- formula(paste(constitColName,"~model(7)"))
  rloadest5param <- loadReg2(
    loadReg(loadRegFormula, data = siteConstitRloadest, 
            flow = qColName, dates = dateColName, time.step = "day",
            flow.units = getInfo(siteMeta, 'flow.units', unit.format = "rloadest"), 
            conc.units = getInfo(siteMeta, 'conc.units', unit.format = "rloadest"),
            load.units = getInfo(siteMeta, 'load.units')), 
    site.id = getInfo(siteMeta, 'site.id'))
  
  interpRect <- loadInterp(
    interp.format = "conc", interp.function = rectangularInterpolation,
    data = siteConstit, metadata = siteMeta)
  
  rloadest5forComp <- loadReg2( # only difference is the data (non-censored)
    loadReg(loadRegFormula, data = siteConstit, 
            flow = qColName, dates = dateColName, time.step = "day",
            flow.units = getInfo(siteMeta, 'flow.units', unit.format = "rloadest"), 
            conc.units = getInfo(siteMeta, 'conc.units', unit.format = "rloadest"),
            load.units = getInfo(siteMeta, 'load.units')), 
    site.id = getInfo(siteMeta, 'site.id'))
  comp <- loadComp(
    reg.model = rloadest5forComp, interp.format = "conc", interp.function = rectangularInterpolation, 
    interp.data = siteConstit)
  
  #list of all model objects
  allModels[[constitSite]] <- list(comp = comp, interpRect = interpRect, 
                                   rloadest5param = rloadest5param)
  
  #make predictions
  pconc_rload <- predictSolute(rloadest5param, "conc", siteQ, se.pred = TRUE, date = TRUE)
  pconc_interp <- predictSolute(interpRect, "conc", siteQ, se.pred = TRUE, date = TRUE)
  pconc_comp <- predictSolute(comp, "conc", siteQ, se.pred = TRUE, date = TRUE)
  pflux_rload <- predictSolute(rloadest5param, "flux", siteQ, se.pred = TRUE, date = TRUE)
  pflux_interp <- predictSolute(interpRect, "flux", siteQ, se.pred = TRUE, date = TRUE)
  pflux_comp <- predictSolute(comp, "flux", siteQ, se.pred = TRUE, date = TRUE)
  nPreds <- nrow(siteQ)
  allPreds <- bind_rows(pconc_rload, pconc_interp, pconc_comp)
  allPreds$model <- rep(c('rloadest','interp','composite'), each=nPreds)
    
  #TODO: model metrics for non-rloadest models
  #combine into DF with row for each model
  #need to extract fitted model so rloadest functions can be used
  metrics <- bind_cols(
    data.frame(summarizeModel(rloadest5param)[1:2]), # site/constit info
    data.frame(REG=summarizeModel(rloadest5param)[-(1:2)]),
    data.frame(INT=summarizeModel(interpRect, irregular.timesteps.ok=TRUE)[-(1:2)]),
    data.frame(CMP=summarizeModel(comp, newdata=siteQ, irregular.timesteps.ok=TRUE)[-(1:2)]))
  write.csv(x = metrics, file = file.path(outputFolder, constitName, "modelMetrics", paste0(constitSite, ".csv")), row.names = FALSE)
  
  #make predictions
  annualSummary <- bind_rows(
    aggregateSolute(pflux_rload, siteMeta, agg.by = "water year", model.name = "rloadest"),
    aggregateSolute(pflux_interp, siteMeta, agg.by = "water year", model.name = "interpolation"),
    aggregateSolute(pflux_comp, siteMeta, agg.by = "water year", model.name = "composite"))
  annualSummary <- reshape(annualSummary, idvar = "water_year", direction = "wide", 
                           v.names = c("Conc","SE", "CI_lower", "CI_upper"), timevar = "model")
  write.csv(x = annualSummary, file = file.path(outputFolder, constitName, "annual", 
                                                paste0(constitSite, '.csv')), row.names=FALSE)
  
  multiYearSummary <- bind_rows(
    aggregateSolute(pflux_rload, siteMeta, agg.by = "mean water year", model.name = "rloadest"),
    aggregateSolute(pflux_interp, siteMeta, agg.by = "mean water year", model.name = "interpolation"),
    aggregateSolute(pflux_comp, siteMeta, agg.by = "mean water year", model.name = "composite"))
  multiYearSummary <- reshape(multiYearSummary, idvar = "Site_Id", direction = "wide", 
                              v.names = c("multi_year_avg", "multiSE", "CI_lower", "CI_upper"), timevar = "model")
  write.csv(x = multiYearSummary, file = file.path(outputFolder, constitName, "multiYear", paste0(constitSite, '.csv')), row.names=FALSE)
  
  #plots
  writePDFreport(file = file.path(outputFolder, constitName, paste(constitSite, "report.pdf", sep = "_")),
                 intdat = siteConstit[1:5], estdat = siteQ, allPreds = allPreds, 
                 meta = siteMeta, inputCSV = inputMetrics, annualCSV = annualPreds)
    
  message(paste('Finished processing constituent file', fileDF$constitFile[i], '\n'))
}

#close the final pdf
dev.off()

allInputs <- summarizeCsvs('inputs', fileDF, outputFolder) 
allAnnual <- summarizeCsvs('annual', fileDF, outputFolder) 
allMultiYear <- summarizeCsvs('multiYear', fileDF, outputFolder) 
allModelMetrics <- summarizeCsvs('modelMetrics', fileDF, outputFolder)

