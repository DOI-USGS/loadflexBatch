# This script runs loadflex in "batch mode": It fits several load estimation 
# models for many constituents at many sites, generates and saves the 
# predictions, and produces summaries over all models, constituents, and sites. 
# See https://github.com/USGS-R/loadflexBatch/blob/master/blog.md for an
# overview and instructions on how to use this file.

#### User inputs ####

inputs <- yaml::yaml.load_file('ANA_sites_170529.yml')


#### Load packages, read inputs, set up directories ####

library(dplyr)
library(loadflex)
library(tools)
library(rloadest)
source('batchHelperFunctions.R') # functions that support this script are stored here

# Read the site info file and attach corresponding files
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
outDetailsDirs <- file.path(rep(outConstitDirs, each=5), rep(c("inputs","annual","multiYear","modelMetrics","plots"), times=nConstits))
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
    
    # Fit the rloadest model (use the smwrQW censoring format)
    # L5: center(log(FLOW)) + center(dectime(DATE)) + fourier(DATE)
    loadRegFormulaL5 <- formula(paste(constitColName,"~model(7)"))
    siteConstitRloadest <- setNames(
      siteConstit[c(dateColName, qwconstitColName, qColName)],
      c(dateColName, constitColName, qColName))
    rloadest5param <- loadReg2(
      loadReg(
        loadRegFormulaL5, data = siteConstitRloadest, 
        flow = qColName, dates = dateColName, time.step = "day",
        flow.units = getUnits(siteMeta, 'flow', format = "rloadest"), 
        conc.units = getUnits(siteMeta, 'conc', format = "rloadest"),
        load.units = getUnits(siteMeta, 'flux', format = "rloadest")), 
      site.id = getInfo(siteMeta, 'site.id'),
      pred.format = 'conc')
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
    
    # Fit the interpolation model (no censoring)
    interpRect <- loadInterp(
      interp.format = "conc", interp.function = rectangularInterpolation,
      data = siteConstit, metadata = siteMeta)
    
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
    
    # Create list of all model objects
    allModels <- list(RL5=rloadest5param, RL7=rloadest7param, CMP=comp, INT=interpRect) #, REGU=rloadest5nocens
    
    
    #### Create output data files ####
    
    # Compute and save info on the site, constituent, and input datasets.
    # Compute num.censored specially here because we're using rloadest format
    # for censored data (smwrQW format) and will eventually have something
    # simpler in place for loadflex, at which point we'll add that to
    # summarizeInputs.
    inputMetrics <- summarizeInputs(siteMeta, fitdat=siteConstit, estdat=siteQ)
    inputMetrics$fitdat.num.censored <- length(which(siteConstit[['status']] == 2))
    inputMetrics$estdat.num.censored <- NULL # assuming no censoring in Q
    inputMetrics <- inputMetrics %>%
      mutate(
        loadflex.version = loadflexVersion, 
        run.date = batchStartTime)
    write.csv(
      x = inputMetrics,
      file = file.path(inputs$outputFolder, constitName, "inputs", paste0(matchingSite, '.csv')),
      row.names=FALSE)
    
    # Summarize each model
    metrics <- bind_cols(
      data.frame(summarizeModel(allModels[[1]])[1:2]), # site.id and constituent columns just once
      lapply(names(allModels), function(mod) { # model-specific columns
        modSum <- switch( 
          class(allModels[[mod]]), # slightly different call for each model type
          loadReg2 = summarizeModel(allModels[[mod]]),
          loadComp = summarizeModel(allModels[[mod]], newdata=siteQ, irregular.timesteps.ok=TRUE),
          loadInterp = summarizeModel(allModels[[mod]], irregular.timesteps.ok=TRUE)
        )
        modSum[-(1:2)] %>% # don't duplicate the site.id and constituent columns
          setNames(paste0(mod, '.', names(.))) # add the "REG.", "CMP.", etc. prefix
      })) %>%
      mutate(
        loadflex.version = loadflexVersion, 
        run.date = batchStartTime)
    write.csv(
      x = metrics, 
      file = file.path(inputs$outputFolder, constitName, "modelMetrics", paste0(matchingSite, ".csv")),
      row.names = FALSE)
    
    # Prepare to convert prediction units if needed
    model.load.rate.units <- getInfo(siteMeta, 'load.rate.units')
    input.load.rate.units <- loadflex:::translateFreeformToUnitted(inputs$loadRateUnits)
    conv.load.rate <- loadflex:::convertUnits(model.load.rate.units, input.load.rate.units)
    
    # Make predictions
    predsLoad <- lapply(allModels, function(mod) {
      (if(is(mod, 'loadComp')) {
        suppressWarnings(predictSolute(mod, "flux", siteQ, se.pred=FALSE, date=TRUE)) %>%
          mutate(se.pred=NA)
      } else {
        predictSolute(mod, "flux", siteQ, se.pred=TRUE, date=TRUE)
      }) %>%
        mutate(
          fit = fit * conv.load.rate,
          se.pred = se.pred * conv.load.rate
        )
    })
    
    # Predict annual fluxes
    message(" * generating annual mean load estimates...")
    annualSummary <- bind_rows(lapply(names(predsLoad), function(mod) {
      message(paste0('   ', mod, '...'), appendLF = FALSE)
      if(is(allModels[[mod]], 'loadReg2')) {
        siteQ %>%
          mutate(Water_Year=smwrBase::waterYear(siteQ[[dateColName]])) %>%
          group_by(Water_Year) %>%
          do({
            # years with a lot of NaNs in the se.pred really slow down the
            # estimation. Remove those rows (changes the )
            siteYearQ <- .
            dailyPreds <- predsLoad[[mod]] %>% 
              mutate(Water_Year=smwrBase::waterYear(predsLoad[[mod]][[dateColName]])) %>%
              filter(Water_Year == siteYearQ$Water_Year[1])
            if(length(which(!is.finite(dailyPreds$se.pred))) > 0) {
              message('skipping NaN-riddled ', as.character(siteYearQ$Water_Year[1]), '...', appendLF = FALSE)
              data.frame(Flux=NaN, SEP=NaN, Ndays=0) # Ndays=0 ensures this year is skipped for total, too
            } else {
              siteYearQ <- filter(siteYearQ, is.finite(dailyPreds$se.pred)) %>%
                select(-Water_Year)
              predLoad(getFittedModel(allModels[[mod]]), newdata=siteYearQ, by='water year', allow.incomplete=TRUE)
            }
          }) %>%
          ungroup() %>%
          mutate(
            Flux_Rate = Flux * conv.load.rate,
            SE = SEP * conv.load.rate,
            n = Ndays,
            CI_lower = L95 * conv.load.rate,
            CI_upper = U95 * conv.load.rate,
            model = mod
          ) %>%
          select(Water_Year, Flux_Rate, SE, n, CI_lower, CI_upper, model)
      } else {
        suppressWarnings(aggregateSolute(
          predsLoad[[mod]], siteMeta, agg.by="water year", format='flux rate')) %>%
          mutate(
            SE = NA,
            CI_lower = NA,
            CI_upper = NA,
            model=mod)
      }
    })) %>%
      mutate(site.id=getInfo(siteMeta, 'site.id'), constituent=getInfo(siteMeta, 'constituent')) %>%
      select(site.id, constituent, model, everything())
    col_order <- c(
      'site.id', 'constituent', 'Water_Year',
      sapply(names(allModels), function(mod) paste0(mod, '.', c('n', 'Flux_Rate', 'SE', 'CI_lower', 'CI_upper'))))
    annualSummary <- annualSummary %>%
      tidyr::gather(var, val, Flux_Rate, SE, n, CI_lower, CI_upper) %>%
      mutate(var=ordered(paste0(model, '.', var))) %>%
      select(-model) %>%
      tidyr::spread(var, val) %>%
      select_(.dots=col_order)
    annualSummary <- annualSummary %>%
      mutate(
        loadflex.version = loadflexVersion, 
        run.date = batchStartTime)
    message('done!')
    write.csv(
      x = annualSummary, 
      file = file.path(inputs$outputFolder, constitName, "annual", paste0(matchingSite, '.csv')),
      row.names=FALSE)
    
    # Predict the multi-year average flux, complete years only
    message(" * generating multi-year mean load estimates...", appendLF = FALSE)
    multiYearSummary <- bind_rows(lapply(names(predsLoad), function(mod) {
      message(paste0(mod, '...'), appendLF=FALSE)
      completeWaterYears <- annualSummary %>% filter(.[[paste0(mod,'.n')]] >= inputs$minDaysPerYear) %>% .$Water_Year
      completeSiteQ <- mutate(siteQ, Water_Year = smwrBase::waterYear(date)) %>%
        filter(Water_Year %in% completeWaterYears)
      (if(is(allModels[[mod]], 'loadReg2')) {
        predLoad(getFittedModel(allModels[[mod]]), newdata=completeSiteQ, by='total', allow.incomplete=TRUE) %>%
          mutate(
            Flux_Rate = Flux * conv.load.rate,
            SE = SEP * conv.load.rate,
            CI_lower = L95 * conv.load.rate,
            CI_upper = U95 * conv.load.rate,
            years.record = length(unique(annualSummary$Water_Year)),            
            years.complete = length(unique(completeWaterYears))
          ) %>%
          select(Flux_Rate, SE, CI_lower, CI_upper, years.record, years.complete)
      } else {
        suppressWarnings(aggregateSolute(
          predsLoad[[mod]], siteMeta, agg.by="mean water year", 
          format='flux rate', min.n=inputs$minDaysPerYear, ci.agg=FALSE, se.agg=FALSE))
      }) %>%
        mutate(model=mod)
    })) %>%
      mutate(site.id=getInfo(siteMeta, 'site.id'), constituent=getInfo(siteMeta, 'constituent')) %>%
      select(site.id, constituent, model, everything())
    multiYearSummary <- reshape(
      multiYearSummary, idvar = c('site.id','constituent','years.record','years.complete'), direction = "wide", 
      v.names = c("Flux_Rate", "SE", "CI_lower", "CI_upper"), timevar = "model")
    multiYearSummary <- setNames( # replace the ".REG" or ".CMP" suffixes with "REG.", "CMP.", prefixes
      multiYearSummary, 
      sub(pattern=sprintf('(.*)\\.(%s)', paste0(names(allModels), collapse='|')), 
          replacement='\\2.\\1', names(multiYearSummary)))
    multiYearSummary <- multiYearSummary %>%
      mutate(
        loadflex.version = loadflexVersion, 
        run.date = batchStartTime)
    message('done!')
    write.csv(
      x = multiYearSummary, 
      file = file.path(inputs$outputFolder, constitName, "multiYear", paste0(matchingSite, '.csv')),
      row.names=FALSE)
    
    #### Create plots  ####
    
    # Start this constituent's pdf file
    graphics.off()
    pdf(height = 10, width = 7.5, paper = "letter",
        file = file.path(inputs$outputFolder, constitName, "plots", sprintf("%s.pdf", siteName)))
    
    # Add plots to the pdf we have open already
    writePDFreport(loadModels = allModels, estdat = siteQ, siteMeta = siteMeta,
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
message("use Adobe or equivalent to combine pdfs")
