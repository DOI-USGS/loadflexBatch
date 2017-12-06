# This script runs loadflex in "batch mode": It fits several load estimation 
# models for many constituents at many sites, generates and saves the 
# predictions, and produces summaries over all models, constituents, and sites. 
# See https://github.com/USGS-R/loadflexBatch/blob/master/blog.md for an
# overview and instructions on how to use this file.

#### User inputs ####

control_file <- 'three_ANA_sites.yml'


#### Load packages & code files, read inputs, set up directories ####

library(dplyr)
library(loadflex)
library(tools)
library(rloadest)
source('batchHelperFunctions.R') # functions that support this script are stored here
bealesFiles <- sapply(file.path('batch_Beales', c('altmod.R','collapse_stratbins.R','mod.R','nsamp.R','predict_ratio.R')), source)

# Read the site info file and attach corresponding files
inputs <- readInputs(control_file)
allSiteInfo <- combineSpecs(inputs)
siteFileSets <- matchFiles(allSiteInfo)

# If outputTimestamp==TRUE, add timestamp to the output folder name
if(inputs$outputTimestamp) {
  inputs$outputFolder <- paste0(inputs$outputFolder, '_', format(Sys.time(), '%y%m%d_%H%M%S'))
}

# Create output directories
constits <- unique(siteFileSets$constituent.CONC)
nConstits <- length(constits)
outConstitDirs <- file.path(rep(inputs$outputFolder, nConstits), constits)
outDetailsDirs <- file.path(rep(outConstitDirs, each=5), rep(c("inputs",inputs$resolutions,"modelMetrics","plots"), times=nConstits))
sapply(outDetailsDirs, dir.create, recursive=TRUE, showWarnings = FALSE)


#### Loop over constituents ####

# Loop over each constituent, creating a pdf of all sites and models for that 
# constituent (plus many smaller, site- and model-specific files)
loadflexVersion <- as.character(packageVersion('loadflex'))
batchStartTime <- Sys.time() # for GitHub repo use "2017-03-22 14:31:59"
message('running loadflex version ', loadflexVersion, ' in batch mode at ', batchStartTime)
for(constitName in constits) { # constitName='NO3'
  
  #### Loop over sites within constituent ####
  
  # Loop over sites having this constituent
  constSites <- filter(siteFileSets, constituent.CONC == constitName)
  for(siteName in constSites$matching.site) { # siteName=constSites$matching.site[2]
    
    # Plan for and announce this iteration
    constitSiteInfo <- filter(constSites, matching.site == siteName)
    matchingSite <- constitSiteInfo$matching.site # this is how we'll name the output files
    message(paste0(
      'processing ', constitName,' at site ', matchingSite,' with files\n',
      ' ', constitSiteInfo$filepath.CONC, ' and\n',
      ' ', constitSiteInfo$filepath.FLOW))
    
    # Read in appropriate files
    siteQ <- read.csv(constitSiteInfo$filepath.FLOW, stringsAsFactors = FALSE)
    siteConstit <- read.csv(constitSiteInfo$filepath.CONC, stringsAsFactors = FALSE)
    if(nrow(siteConstit) == 0) {
      message("empty table at ", constitSiteInfo$filepath.CONC)
      next
    }
    
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
    
    # Convert date text to Dates (inputs must be >= daily for this batch script)
    siteQ[[dateColName]] <- as.Date(siteQ[[dateColName]], format='%Y-%m-%d')
    siteConstit[[dateColName]] <- as.Date(siteConstit[[dateColName]], format='%Y-%m-%d')
    
    # Filter to date ranges from siteInfo.csv if specified
    if(!is.na(constitSiteInfo$date.start.FLOW)) siteQ <- siteQ[which(siteQ[[dateColName]] >= constitSiteInfo$date.start.FLOW),]
    if(!is.na(constitSiteInfo$date.end.FLOW)) siteQ <- siteQ[which(siteQ[[dateColName]] <= constitSiteInfo$date.end.FLOW),]
    if(!is.na(constitSiteInfo$date.start.CONC)) siteConstit <- siteConstit[which(siteConstit[[dateColName]] >= constitSiteInfo$date.start.CONC),]
    if(!is.na(constitSiteInfo$date.end.CONC)) siteConstit <- siteConstit[which(siteConstit[[dateColName]] <= constitSiteInfo$date.end.CONC),]
    
    # Deal with different discharge/consituent drainage areas. Compute discharge
    # for both the Q and constit data.frames as if at constituent site
    if(constitSiteInfo$basin.area.FLOW != constitSiteInfo$basin.area.CONC) {
      QratioCF <- constitSiteInfo$basin.area.CONC/constitSiteInfo$basin.area.FLOW
      siteConstit[[qColName]] <- siteConstit[[qColName]]*QratioCF
      siteQ[[qColName]] <- siteQ[[qColName]]*QratioCF
      message(sprintf(" * scaling discharge by basin area: multiplying by %1.3f", QratioCF))
    }
    
    # Remove duplicate observations
    constitDupes <- table(siteConstit[[dateColName]]) %>% .[.>1] %>% names()
    nonFirstDupes <- unlist(lapply(constitDupes, function(cD) {
      dupes <- which(as.character(siteConstit[[dateColName]]) == cD)
      return(dupes[-1])
    }))
    if(length(nonFirstDupes) > 0) {
      message("  * removing ", length(nonFirstDupes), " duplicate dates (", paste0(constitDupes, collapse=', '), ")")
      siteConstit <- siteConstit[-nonFirstDupes,]
    }
    
    # Format censored data for rloadest. For ANA, "Status 0 means null or blank. 
    # Status 1 means a valid value, and 2 means that the respective value is a 
    # detection limit."
    censor.statuses <- c('0'='', '1'='', '2'='<')
    siteConstit[[qwconstitColName]] <- 
      smwrQW::as.lcens(
        values=ifelse(siteConstit[['status']] == 0, NA, siteConstit[[constitColName]]), 
        detlim=ifelse(siteConstit[['status']] < 2, 0, siteConstit[[constitColName]]), 
        censor.codes=censor.statuses[as.character(siteConstit[['status']])])
    siteConstit[[constitColName]] <- 
      ifelse(
        siteConstit[['status']] == 0, 
        NA, 
        ifelse(
          siteConstit[['status']] == 2, 
          siteConstit[[constitColName]] / 2, # it's still bad, but use half MDL because better than full MDL
          siteConstit[[constitColName]]))
    # Remove NA values, some of which may have been added by checking status
    siteConstit <- siteConstit[!is.na(siteConstit[[constitColName]]), ]
    
    #### Create loadflex models ####
    
    # Create a formal metadata object
    siteMeta <- metadata(
      constituent = constitColName, consti.name = constitColName, conc.units = constitSiteInfo$units.CONC, 
      flow = qColName, flow.units = constitSiteInfo$units.FLOW, load.units = inputs$loadUnits, 
      load.rate.units = paste(inputs$loadUnits, 'd^-1'),  # we'll use inputs$loadRateUnits at prediction time
      dates = dateColName,
      site.name = constitSiteInfo$site.name.CONC, site.id = constitSiteInfo$site.id.CONC, 
      lat = constitSiteInfo$lat.CONC, lon = constitSiteInfo$lon.CONC, basin.area = constitSiteInfo$basin.area.CONC,
      flow.site.name = constitSiteInfo$site.name.FLOW, flow.site.id = constitSiteInfo$site.id.FLOW, 
      flow.lat = constitSiteInfo$lat.FLOW, flow.lon = constitSiteInfo$lon.FLOW, flow.basin.area = constitSiteInfo$basin.area.FLOW
    )
    
    # Fix the random number generator seed so we get the same results each time we
    # run the script
    set.seed(9451)
    
    # Create list of all model objects
    allModels <- setNames(as.list(rep(NA, length(inputs$models))), inputs$models)
    # list(RL5=rloadest5param, RL7=rloadest7param, CMP=comp, INT=interpRect, BRE=beales)[inputs$models] #, REGU=rloadest5nocens
    
    if(any(c('RL5','RL7') %in% inputs$models)) {
      # Fit the rloadest model[s] (use the smwrQW censoring format)
      siteConstitRloadest <- setNames(
        siteConstit[c(dateColName, qwconstitColName, qColName)],
        c(dateColName, constitColName, qColName))
      if('RL5' %in% inputs$models) {
        # L5: center(log(FLOW)) + center(dectime(DATE)) + fourier(DATE)
        loadRegFormulaL5 <- formula(paste(constitColName,"~model(7)"))
        rloadest5param <- loadReg2(
          loadReg(
            loadRegFormulaL5, data = siteConstitRloadest, 
            flow = qColName, dates = dateColName, time.step = "day",
            flow.units = getUnits(siteMeta, 'flow', format = "rloadest"), 
            conc.units = getUnits(siteMeta, 'conc', format = "rloadest"),
            load.units = getUnits(siteMeta, 'flux', format = "rloadest")), 
          site.id = getInfo(siteMeta, 'site.id'),
          pred.format = 'conc')
        allModels[['RL5']] <- rloadest5param
      }
      
      if('RL7' %in% inputs$models) {
        # L7: quadratic(log(FLOW)) + quadratic(dectime(DATE)) + fourier(DATE)
        loadRegFormulaL7 <- formula(paste(constitColName,"~model(9)"))
        rloadest7param <- loadReg2(
          loadReg(
            loadRegFormulaL7, data = siteConstitRloadest, 
            flow = qColName, dates = dateColName, time.step = "day",
            flow.units = getUnits(siteMeta, 'flow', format = "rloadest"), 
            conc.units = getUnits(siteMeta, 'conc', format = "rloadest"),
            load.units = getUnits(siteMeta, 'flux', format = "rloadest")), 
          site.id = getInfo(siteMeta, 'site.id'),
          pred.format = 'conc')
        allModels[['RL7']] <- rloadest7param
      }
    }
    
    if('INT' %in% inputs$models) {
      # Fit the interpolation model (no censoring)
      interpRect <- loadInterp(
        interp.format = "conc", interp.function = rectangularInterpolation,
        data = siteConstit, metadata = siteMeta)
      allModels[['INT']] <- interpRect
    }
    
    if('CMP' %in% inputs$models) {
      # Fit the composite model (with rloadest model that doesn't do censoring)
      rloadest5nocens <- loadReg2( # only difference is the data (non-censored)
        loadReg(
          loadRegFormulaL5, data = siteConstit, 
          flow = qColName, dates = dateColName, time.step = "day",
          flow.units = getUnits(siteMeta, 'flow', format = "rloadest"), 
          conc.units = getUnits(siteMeta, 'conc', format = "rloadest"),
          load.units = getUnits(siteMeta, 'flux', format = "rloadest")), 
        site.id = getInfo(siteMeta, 'site.id'),
        pred.format = 'conc')
      comp <- loadComp(
        reg.model = rloadest5nocens, interp.format = "conc", interp.function = rectangularInterpolation, 
        interp.data = siteConstit, store=c('data','fitting.function')) # leave out store='uncertainty' to save 30 secs
      allModels[['CMP']] <- comp
    }
    
    if('BRE' %in% inputs$models) {
      # Check for units assumptions of Beale's ratio estimator implementation
      if(siteMeta@conc.units != 'mg L^-1') stop("For Beale's ratio estimator, constituent units (in siteInfo file) must be 'mg L^-1'")
      if(siteMeta@flow.units != 'm^3 s^-1') stop("For Beale's ratio estimator, flow units (in siteInfo file) must be m^3 s^-1'")
      if(inputs$loadUnits != 'kg') stop("For Beale's ratio estimator, loadUnits (in .yml) must be 'kg'")
      if(loadflex:::translateFreeformToUnitted(inputs$loadRateUnits) != 'kg y^-1') stop("For Beale's ratio estimator, loadRateUnits (in .yml) must be 'kg y^-1'")
      
      # Beale's ratio doesn't have a model object; here we're doing the fitting,
      # diagnostics, and multi-year prediction all at once. The estimator code always generates predictions in kg/y
      beales <- predict_ratio(
        siteQ, siteConstit, # data
        inputs$minDaysPerYear,
        waterYear=TRUE,
        dateName=dateColName,
        qName=qColName,
        constitName=constitColName,
        hi_flow_percentile=80, #Default threshold for designating high-flow observations,
        ratio_strata_nsamp_threshold=10, #Default minimum number of observations required for inclusion of a stratum in the ratio estimate,
        concTrans=1, #constant transformation factor for converting to units of mg/L. NA=1 implies input is mg/L
        qTrans=35.31466) #constant transformation factor for converting Q to units of ft3/s. 35.3 implies input is cms
      class(beales) <- 'loadBeale'
      #   rload_NO3_kg/y serload_NO3_kg/y nstrata
      # 1       88260.66         7910.503       4
      allModels[['BRE']] <- beales
    }
    
    
    #### Create output data files ####
    
    # Summarize the input data
    inputMetrics <- summarizeBatchInputs(siteMeta, siteConstit, siteQ, loadflexVersion, batchStartTime)
    write.csv(
      x = inputMetrics,
      file = file.path(inputs$outputFolder, constitName, "inputs", paste0(matchingSite, '.csv')),
      row.names=FALSE)
    
    # Summarize each model
    metrics <- summarizeMetrics(allModels, siteMeta, loadflexVersion, batchStartTime)
    write.csv(
      x = metrics, 
      file = file.path(inputs$outputFolder, constitName, "modelMetrics", paste0(matchingSite, ".csv")),
      row.names = FALSE)
    
    # Prepare to convert prediction units if needed
    model.load.rate.units <- getInfo(siteMeta, 'load.rate.units')
    input.load.rate.units <- loadflex:::translateFreeformToUnitted(inputs$loadRateUnits)
    conv.load.rate <- loadflex:::convertUnits(model.load.rate.units, input.load.rate.units)
    
    # Predict daily fluxes
    predsLoad <- summarizeDaily(allModels, siteQ, conv.load.rate)
    
    # Predict monthly fluxes
    if('monthly' %in% inputs$resolutions) {
      monthlySummary <- summarizeMonthly(allModels, predsLoad, inputs, siteQ, conv.load.rate, loadflexVersion, batchStartTime)
      write.csv(
        x = monthlySummary, 
        file = file.path(inputs$outputFolder, constitName, "monthly", paste0(matchingSite, '.csv')),
        row.names=FALSE)
    }
    
    # Predict annual fluxes
    if(any(c('annual','multiYear') %in% inputs$resolutions)) {
      annualSummary <- summarizeAnnual(allModels, predsLoad, inputs, siteQ, conv.load.rate, loadflexVersion, batchStartTime)
      if('annual' %in% inputs$resolutions) {
        write.csv(
          x = annualSummary, 
          file = file.path(inputs$outputFolder, constitName, "annual", paste0(matchingSite, '.csv')),
          row.names=FALSE)
      }
    }
    
    # Predict the multi-year average flux, complete years only
    if('multiYear' %in% inputs$resolutions) {
      multiYearSummary <- summarizeMultiYear(allModels, predsLoad, annualSummary, inputs, siteQ, conv.load.rate, loadflexVersion, batchStartTime)
      write.csv(
        x = multiYearSummary, 
        file = file.path(inputs$outputFolder, constitName, "multiYear", paste0(matchingSite, '.csv')),
        row.names=FALSE)
    }    
    
    #### Create plots  ####
    
    # Start this constituent's pdf file
    graphics.off()
    pdf(height = 10, width = 7.5, paper = "letter",
        file = file.path(inputs$outputFolder, constitName, "plots", sprintf("%s.pdf", siteName)))
    
    # Add plots to the pdf we have open already
    writePDFreport(loadModels = allModels, fitdat = siteConstit, estdat = siteQ, siteMeta = siteMeta,
                   loadflexVersion = loadflexVersion, batchStartTime = batchStartTime)
    
    # Close this constituent's pdf file
    dev.off()
    
  }
}

#### Combine outputs from all sites ####

# Combine the outputs from each site-constituent combination into a single table
# per output type and constituent
allInputs <- summarizeCsvs('inputs', siteFileSets, inputs$outputFolder) 
allAnnual <- summarizeCsvs('annual', siteFileSets, inputs$outputFolder) 
allMultiYear <- summarizeCsvs('multiYear', siteFileSets, inputs$outputFolder) 
allModelMetrics <- summarizeCsvs('modelMetrics', siteFileSets, inputs$outputFolder)
# summarizePlots(siteFileSets, inputs$outputFolder) # doesn't work on my computer
message("use Adobe or equivalent to combine pdfs into CONST_plots.pdf")
