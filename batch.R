# This script runs loadflex in "batch mode": It fits several load estimation 
# models for many constituents at many sites, generates and saves the 
# predictions, and produces summaries over all models, constituents, and sites. 
# See https://github.com/USGS-R/loadflexBatch/blob/master/blog.md for an
# overview and instructions on how to use this file.

#------------------User Inputs--------------------#

inputs <- yaml::yaml.load_file('Hirsch_sites.yml')


#-------------------------Load packages, check files, set up directories-----------------------# 

library(dplyr)
library(loadflex)
library(tools)
library(rloadest)
source('batchHelperFunctions.R') # functions that support this script are stored here

# Require at least rloadest 0.4.4 for formula fix
if(compareVersion(as.character(packageVersion("rloadest")),"0.4.4") == -1) {
  stop("rloadest version 0.4.4 or greater is required")
}

# Read the directory structure and site info file to determine what's available
fileDF <- makeFileDF(inputs$inputFolder, constits = inputs$constituents, discharge.folder = inputs$dischargeFolder)
allSiteInfo <- read.csv(file.path(inputs$inputFolder, inputs$siteInfo), stringsAsFactors = FALSE)

# Create output directories
nConstits <- length(inputs$constituents)
outConstitDirs <- file.path(rep(inputs$outputFolder, nConstits), inputs$constituents)
outDetailsDirs <- file.path(rep(outConstitDirs, each=4), rep(c("inputs","annual","multiYear","modelMetrics"), times=nConstits))
sapply(outDetailsDirs, dir.create, recursive=TRUE, showWarnings = FALSE)

lastConstit <- NULL
graphics.off() # we don't want open PDF connections

#-----------------loadflex--------------#

# Loop over unique site-constituent combinations, creating a set of output files
# for each
for(i in 1:nrow(fileDF)) {
  message(paste0('processing constituent file ', fileDF$constitFile[i]))
  
  #TODO: if a constituent file is missing, skip it and keep going (#142)
  
  # Read in appropriate files
  siteQ <- read.csv(fileDF$qFile[i], stringsAsFactors = FALSE)
  siteConstit <- read.csv(fileDF$constitFile[i], stringsAsFactors = FALSE)
  
  # Convert date text to Dates 
  siteConstit$date <- as.Date(siteConstit$date, format='%Y-%m-%d')
  siteQ$date <- as.Date(siteQ$date, format='%Y-%m-%d')
  
  # Pull out appropriate rows of allSiteInfo for Q and constit. Need to extract 
  # constit and sites from file paths so we know what row of site info to check
  constitSite <- basename(file_path_sans_ext(fileDF$constitFile[i])) 
  constitName <- basename(dirname(fileDF$constitFile[i]))
 
  # If switching to a new consituent, open a new pdf. It's important that fileDF
  # is sorted by consituent! We've done this by structuring makeFileDF so it
  # loops over folders
  if(is.null(lastConstit)) {
    pdf(height = 11, width = 8.5, 
        file = file.path(inputs$outputFolder, constitName, sprintf("%s_plots.pdf", constitName)))
    lastConstit <- constitName
  }
  if(constitName != lastConstit) {
    lastConstit <- constitName
    dev.off()
    pdf(height = 11, width = 8.5, 
        file = file.path(inputs$outputFolder, constitName, sprintf("%s_plots.pdf", constitName)))
  }
  
  # Isolate the site, constituent, and flow metadata relevant to this iteration
  constitSiteInfo <- filter(allSiteInfo, matching.site == constitSite, constituent == constitName)
  qSiteInfo <- filter(allSiteInfo, matching.site == constitSite, constituent == 'Q')
  
  # Deal with different discharge/consituent drainage areas
  if(qSiteInfo$basin.area != constitSiteInfo$basin.area) {
    ratio <- constitSiteInfo$basin.area/qSiteInfo$basin.area
    # modify discharge for both DFs
    siteConstit$Q <- siteConstit$Q*ratio
    siteQ$Q <- siteQ$Q*ratio
    message(sprintf(" * scaling discharge by basin area: multiplying by %1.3f", ratio))
  }
  
  # Create metadata. We're depending on consistent column positions to get names
  # and could probably do better (see #190)
  constitColName <- names(siteConstit)[3]
  qwconstitColName <- paste0(constitColName, '_qw')
  qColName <- names(siteConstit)[2]
  dateColName <- names(siteConstit)[1]
  
  # Remove duplicate observations
  constitDupes <- table(siteConstit$date) %>% .[.>1] %>% names()
  nonFirstDupes <- unlist(lapply(constitDupes, function(cD) {
    dupes <- which(as.character(siteConstit$date) == cD)
    return(dupes[-1])
  }))
  if(length(nonFirstDupes) > 0) siteConstit <- siteConstit[-nonFirstDupes,]
  
  # Format censored data for rloadest. For ANA, Status 0 means null or blank. 
  # Status 1 means a valid value, and 2 means that the respective value is a 
  # detection limit."
  censor.statuses <- c('0'='', '1'='', '2'='<')
  siteConstit[[qwconstitColName]] <- 
    smwrQW::as.lcens(
      values=ifelse(siteConstit[['status']] == 0, NA, siteConstit[[constitColName]])/2, 
      detlim=ifelse(siteConstit[['status']] < 2, 0, siteConstit[[constitColName]]), 
      censor.codes=censor.statuses[as.character(siteConstit[['status']])])
  siteConstit[[constitColName]] <- ifelse(
    siteConstit[['status']] == 0, NA, 
    ifelse(siteConstit[['status']] == 2, siteConstit[[constitColName]] / 2, # it's still bad, but use half MDL because better than full MDL
           siteConstit[[constitColName]]))
  # Remove NA values, some of which may have been added by checking status
  siteConstit <- siteConstit[!is.na(siteConstit[[constitColName]]), ]
  
  # Create a formal metadata object. site.id and flow.site.id must both equal
  # constitSite for our input file scheme to work
  siteMeta <- metadata(
    constituent = constitColName, consti.name = constitColName, conc.units = constitSiteInfo$units, 
    flow = qColName, flow.units = qSiteInfo$units, 
    load.units = inputs$loadUnits, load.rate.units = inputs$loadRateUnits, dates = dateColName,
    site.name = constitSiteInfo$site.name, site.id = constitSiteInfo$site.id, lat = constitSiteInfo$lat, lon = constitSiteInfo$lon, basin.area = constitSiteInfo$basin.area,
    flow.site.name = qSiteInfo$site.name, flow.site.id = qSiteInfo$site.id, flow.lat = qSiteInfo$lat, flow.lon = qSiteInfo$lon, flow.basin.area = qSiteInfo$basin.area
  )
  
  # Compute and save info on the site, constituent, and input datasets (we'll 
  # recombine in the next loop). Compute num.censored specially here because 
  # we're using rloadest format for censored data (smwrQW format) and will 
  # eventually have something simpler in place for loadflex, at which point
  # we'll add that to summarizeInputs.
  inputMetrics <- summarizeInputs(siteMeta, fitdat=siteConstit, estdat=siteQ)
  inputMetrics$fitdat.num.censored <- length(which(!is.na(siteConstit[[qwconstitColName]]@.Data[,'detlim'])))
  inputMetrics$estdat.num.censored <- NULL # assuming there isn't and shouldn't be censoring in Q. is that right?
  write.csv(inputMetrics, file.path(inputs$outputFolder, constitName, "inputs", paste0(constitSite, '.csv')), row.names=FALSE)
  
  # Fix the random number generator seed so we get the same results each time we
  # run the script
  set.seed(9451)
  
  # Fit the rloadest model
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
  
  # Fit the interpolation model (no censoring)
  interpRect <- loadInterp(
    interp.format = "conc", interp.function = rectangularInterpolation,
    data = siteConstit, metadata = siteMeta)
  
  # Fit the composite model (with rloadest model that doesn't do censoring)
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
  
  # Create list of all model objects
  allModels <- list(composite=comp, interp=interpRect, rloadest=rloadest5param)
  
  # Make predictions
  pconc_rload <- predictSolute(rloadest5param, "conc", siteQ, se.pred = TRUE, date = TRUE)
  pconc_interp <- predictSolute(interpRect, "conc", siteQ, se.pred = TRUE, date = TRUE)
  pconc_comp <- predictSolute(comp, "conc", siteQ, se.pred = TRUE, date = TRUE)
  pflux_rload <- predictSolute(rloadest5param, "flux", siteQ, se.pred = TRUE, date = TRUE)
  pflux_interp <- predictSolute(interpRect, "flux", siteQ, se.pred = TRUE, date = TRUE)
  pflux_comp <- predictSolute(comp, "flux", siteQ, se.pred = TRUE, date = TRUE)
  nPreds <- nrow(siteQ)
  allConcPreds <- bind_rows(pconc_rload, pconc_interp, pconc_comp)
  allConcPreds$model <- rep(c('rloadest','interp','composite'), each=nPreds)
    
  # Summarize each model
  metrics <- bind_cols(
    data.frame(summarizeModel(rloadest5param)[1:2]), # site/constit info
    data.frame(REG=summarizeModel(rloadest5param)[-(1:2)]),
    data.frame(INT=summarizeModel(interpRect, irregular.timesteps.ok=TRUE)[-(1:2)]),
    data.frame(CMP=summarizeModel(comp, newdata=siteQ, irregular.timesteps.ok=TRUE)[-(1:2)]))
  write.csv(x = metrics, file = file.path(inputs$outputFolder, constitName, "modelMetrics", paste0(constitSite, ".csv")), row.names = FALSE)
  
  # Predict annual fluxes
  annualSummary <- bind_rows(
    mutate(aggregateSolute(pflux_rload, siteMeta, agg.by = "water year"), model = "REG"),
    mutate(aggregateSolute(pflux_interp, siteMeta, agg.by = "water year"), model = "INT"),
    mutate(aggregateSolute(pflux_comp, siteMeta, agg.by = "water year"), model = "CMP")) %>%
    mutate(site.id=siteMeta@site.id, constituent=siteMeta@constituent) %>%
    select(site.id, constituent, model, everything())
  annualSummary <- reshape(
    annualSummary, idvar = "water_year", direction = "wide", 
    v.names = c("Conc","SE", "CI_lower", "CI_upper"), timevar = "model")
  annualSummary <- setNames(
    annualSummary, 
    sub(pattern='(.*)\\.(REG|INT|CMP)', replacement='\\2.\\1', names(annualSummary)))
  write.csv(
    x = annualSummary, 
    file = file.path(inputs$outputFolder, constitName, "annual", paste0(constitSite, '.csv')), row.names=FALSE)
  
  # Predict the multi-year average flux
  multiYearSummary <- bind_rows(
    mutate(aggregateSolute(pflux_rload, siteMeta, agg.by = "mean water year"), model = "REG"),
    mutate(aggregateSolute(pflux_interp, siteMeta, agg.by = "mean water year"), model = "INT"),
    mutate(aggregateSolute(pflux_comp, siteMeta, agg.by = "mean water year"), model = "CMP")) %>%
    mutate(site.id=siteMeta@site.id, constituent=siteMeta@constituent) %>%
    select(site.id, constituent, model, everything())
  multiYearSummary <- reshape(
    multiYearSummary, direction = "wide", idvar=c("site.id", "constituent"),
    v.names = c("Multi_Year_Avg", "multiSE", "CI_lower", "CI_upper", "years.record", "years.complete"), timevar = "model")
  multiYearSummary <- setNames(
    multiYearSummary, 
    sub(pattern='(.*)\\.(REG|INT|CMP)', replacement='\\2.\\1', names(multiYearSummary)))
  write.csv(
    x = multiYearSummary, 
    file = file.path(inputs$outputFolder, constitName, "multiYear", paste0(constitSite, '.csv')), row.names=FALSE)
  
  # Add plots to the pdf we have open for writing
  writePDFreport(
    file = file.path(inputs$outputFolder, constitName, paste(constitSite, "report.pdf", sep = "_")),
    load.models = allModels, estdat = siteQ, siteMeta = siteMeta)
}

# Close the final pdf
dev.off()

# Combine the outputs from each site-constituent combination into a single table
# per output type and constituent
allInputs <- summarizeCsvs('inputs', fileDF, inputs$outputFolder) 
allAnnual <- summarizeCsvs('annual', fileDF, inputs$outputFolder) 
allMultiYear <- summarizeCsvs('multiYear', fileDF, inputs$outputFolder) 
allModelMetrics <- summarizeCsvs('modelMetrics', fileDF, inputs$outputFolder)
