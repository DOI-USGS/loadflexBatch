# This script runs loadflex in "batch mode": It fits several load estimation 
# models for many constituents at many sites, generates and saves the 
# predictions, and produces summaries over all models, constituents, and sites. 
# See https://github.com/USGS-R/loadflexBatch/blob/master/blog.md for an
# overview and instructions on how to use this file.

#### User inputs ####

inputs <- yaml::yaml.load_file('Hirsch_sites.yml')


#### Load packages, read inputs, set up directories ####

library(dplyr)
library(loadflex)
library(tools)
library(rloadest)
source('batchHelperFunctions.R') # functions that support this script are stored here

# Read the site info file and attach corresponding files
allSiteInfo <- combineSpecs(inputs)
siteFileSets <- matchFiles(allSiteInfo)

# Create output directories
constits <- unique(siteFileSets$constituent.CONC)
nConstits <- length(constits)
outConstitDirs <- file.path(rep(inputs$outputFolder, nConstits), constits)
outDetailsDirs <- file.path(rep(outConstitDirs, each=4), rep(c("inputs","annual","multiYear","modelMetrics"), times=nConstits))
sapply(outDetailsDirs, dir.create, recursive=TRUE, showWarnings = FALSE)


#### Loop over constituents ####

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
  
  
  # Remove duplicate observations
  constitDupes <- table(siteConstit$date) %>% .[.>1] %>% names()
  nonFirstDupes <- unlist(lapply(constitDupes, function(cD) {
    dupes <- which(as.character(siteConstit$date) == cD)
    return(dupes[-1])
  }))
  if(length(nonFirstDupes) > 0) {
    message("  * removing ", length(nonFirstDupes), " rows with duplicate dates")
    siteConstit <- siteConstit[-nonFirstDupes,]
# Loop over each constituent, creating a pdf of all sites and models for that 
# constituent (plus many smaller, site- and model-specific files)
for(constitName in constits) {
  
  # Start this constituent's pdf file
  graphics.off()
  pdf(height = 11, width = 8.5, 
      file = file.path(inputs$outputFolder, constitName, sprintf("%s_plots.pdf", constitName)))
  
  #### Loop over sites within constituent ####
  
  # Loop over sites having this constituent
  constSites <- filter(siteFileSets, constituent.CONC == constitName)
  for(siteName in constSites$matching.site) {
    # Identify and require column names as given in site info
    constitColName <- constitSiteInfo$constituent.CONC
    qwconstitColName <- paste0(constitColName, '_qw')
    qColName <- constitSiteInfo$constituent.FLOW
    dateColName <- inputs$date
    missingConstitCols <- names(which(sapply(c(dateColName, constitColName, qColName), function(col) !(col %in% colnames(siteConstit)))))
    missingQCols <- names(which(sapply(c(dateColName, qColName), function(col) !(col %in% colnames(siteQ)))))
    if(length(missingConstitCols) > 0) stop("missing these columns in the constituent data: ", paste0(missingConstitCols, collapse=', '))
    if(length(missingQCols) > 0) stop("missing these columns in the discharge data: ", paste0(missingQCols, collapse=', '))
    #### Munge the input data ####
    #### Create loadflex models ####
    #### Create output data files ####
    #### Create plots  ####
  }
  
  # Format censored data for rloadest. For ANA, Status 0 means null or blank. 
  # Status 1 means a valid value, and 2 means that the respective value is a 
  # detection limit."
  censor.statuses <- c('0'='', '1'='', '2'='<')
  siteConstit[[qwconstitColName]] <- 
    smwrQW::as.lcens(
      values=ifelse(siteConstit[['status']] == 0, NA, siteConstit[[constitColName]]), 
      detlim=ifelse(siteConstit[['status']] < 2, 0, siteConstit[[constitColName]]), 
      censor.codes=censor.statuses[as.character(siteConstit[['status']])])
  siteConstit[[constitColName]] <- 
    ifelse(siteConstit[['status']] == 0, 
           NA, 
           ifelse(siteConstit[['status']] == 2, 
                  siteConstit[[constitColName]] / 2, # it's still bad, but use half MDL because better than full MDL
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
            flow.units = getUnits(siteMeta, 'flow', format = "rloadest"), 
            conc.units = getUnits(siteMeta, 'conc', format = "rloadest"),
            load.units = getUnits(siteMeta, 'flux', format = "rloadest")), 
    site.id = getInfo(siteMeta, 'site.id'),
    pred.format = 'conc')
  
  # Fit the interpolation model (no censoring)
  interpRect <- loadInterp(
    interp.format = "conc", interp.function = rectangularInterpolation,
    data = siteConstit, metadata = siteMeta)
  
  # Fit the composite model (with rloadest model that doesn't do censoring)
  rloadest5forComp <- loadReg2( # only difference is the data (non-censored)
    loadReg(loadRegFormula, data = siteConstit, 
            flow = qColName, dates = dateColName, time.step = "day",
            flow.units = getUnits(siteMeta, 'flow', format = "rloadest"), 
            conc.units = getUnits(siteMeta, 'conc', format = "rloadest"),
            load.units = getUnits(siteMeta, 'flux', format = "rloadest")), 
    site.id = getInfo(siteMeta, 'site.id'))
  comp <- loadComp(
    reg.model = rloadest5forComp, interp.format = "conc", interp.function = rectangularInterpolation, 
    interp.data = siteConstit)
  
  # Create list of all model objects
  allModels <- list(CMP=comp, INT=interpRect, REGC=rloadest5param, REGU=rloadest5forComp)
  
  # Make predictions
  predsLoad <- lapply(allModels, predictSolute, "flux", siteQ, se.pred = TRUE, date = TRUE)
  predsConc <- lapply(allModels, predictSolute, "conc", siteQ, se.pred = TRUE, date = TRUE)
    
  # Summarize each model
  metrics <- bind_cols(
    data.frame(summarizeModel(rloadest5param)[1:2]), # site/constit info
    data.frame(REGC=summarizeModel(rloadest5param)[-(1:2)]),
    data.frame(REGU=summarizeModel(rloadest5forComp)[-(1:2)]),
    data.frame(INT=summarizeModel(interpRect, irregular.timesteps.ok=TRUE)[-(1:2)]),
    data.frame(CMP=summarizeModel(comp, newdata=siteQ, irregular.timesteps.ok=TRUE)[-(1:2)]))
  write.csv(x = metrics, file = file.path(inputs$outputFolder, constitName, "modelMetrics", paste0(constitSite, ".csv")), row.names = FALSE)
  
  # Predict annual fluxes
  annualSummary <- bind_rows(lapply(names(predsLoad), function(mod) {
    preds <- predsLoad[[mod]]
    mutate(aggregateSolute(preds, siteMeta, agg.by="water year", format='flux rate'), model=mod)
  })) %>%
    mutate(site.id=siteMeta@site.id, constituent=siteMeta@constituent) %>%
    select(site.id, constituent, model, everything())
  annualSummary <- reshape(
    annualSummary, idvar = c('site.id','constituent',"Water_Year"), direction = "wide", 
    v.names = c("Flux_Rate", "SE", "CI_lower", "CI_upper"), timevar = "model")
  annualSummary <- setNames(
    annualSummary, 
    sub(pattern='(.*)\\.(REGC|REGU|INT|CMP)', replacement='\\2.\\1', names(annualSummary)))
  write.csv(
    x = annualSummary, 
    file = file.path(inputs$outputFolder, constitName, "annual", paste0(constitSite, '.csv')), row.names=FALSE)
  
  # Predict the multi-year average flux
  multiYearSummary <- bind_rows(lapply(names(predsLoad), function(mod) {
    preds <- predsLoad[[mod]]
    mutate(aggregateSolute(preds, siteMeta, agg.by="mean water year", format='flux rate'), model=mod)
  })) %>%
    mutate(site.id=siteMeta@site.id, constituent=siteMeta@constituent) %>%
    select(site.id, constituent, model, everything())
  multiYearSummary <- reshape(
    multiYearSummary, idvar = c('site.id','constituent'), direction = "wide", 
    v.names = c("Flux_Rate", "SE", "CI_lower", "CI_upper", "years.record", "years.complete"), timevar = "model")
  multiYearSummary <- setNames(
    multiYearSummary, 
    sub(pattern='(.*)\\.(REGC|REGU|INT|CMP)', replacement='\\2.\\1', names(annualSummary)))
  write.csv(
    x = multiYearSummary, 
    file = file.path(inputs$outputFolder, constitName, "multiYear", paste0(constitSite, '.csv')), row.names=FALSE)
  
  # Add plots to the pdf we have open for writing
  writePDFreport(
    file = file.path(inputs$outputFolder, constitName, paste(constitSite, "report.pdf", sep = "_")),
    load.models = allModels, estdat = siteQ, siteMeta = siteMeta)
  # Close this constituent's pdf file
  dev.off()
}

# Close the final pdf
dev.off()
#### Combine outputs from all sites ####

# Combine the outputs from each site-constituent combination into a single table
# per output type and constituent
allInputs <- summarizeCsvs('inputs', fileDF, inputs$outputFolder) 
allAnnual <- summarizeCsvs('annual', fileDF, inputs$outputFolder) 
allMultiYear <- summarizeCsvs('multiYear', fileDF, inputs$outputFolder) 
allModelMetrics <- summarizeCsvs('modelMetrics', fileDF, inputs$outputFolder)
