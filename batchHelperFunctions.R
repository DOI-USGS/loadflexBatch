#' Read and validate the YAML control file
#' @param control_file
readInputs <- function(control_file) {
  inputs <- yaml::yaml.load_file(control_file)
  
  expected <- c(
    "inputFolder","constituents","discharge","date","siteInfo",
    "models","resolutions","minDaysPerYear","regMaxNaNsPerMonth","regMaxNaNsPerSeason","regMaxNaNsPerYear",
    "regBaseYear",
    "loadUnits","loadRateUnits",
    "outputFolder","outputTimestamp")
  miss <- setdiff(expected, names(inputs))
  if(length(miss) > 0) stop(paste0("missing fields in control file '", control_file,"': ", paste0("'", miss, "'", collapse=", ")))
  extra <- setdiff(names(inputs), expected)
  if(length(extra) > 0) stop(paste0("missing fields in control file '", control_file,"': ", paste0("'", extra, "'", collapse=", ")))
  
  return(inputs)
}

#' Make a data.frame describing the constituent and flow sites and the data
#' files for those variables
#' 
#' @param inputs list of user inputs as given in a yml file
combineSpecs <- function(inputs) {
  
  # read in the siteInfo file
  siteInfo <- read.csv(file.path(inputs$inputFolder, inputs$siteInfo), stringsAsFactors = FALSE)
  
  # check for the expected columns
  expectedcols <- c("matching.site", "site.id", "site.name", "lat", "lon", "basin.area", "constituent", "consti.name", "units", "date.start", "date.end")
  misscol <- setdiff(expectedcols, names(siteInfo))
  if(length(misscol) > 0) stop(paste("missing columns in siteInfo file:", paste0("'", misscol, "'", collapse=", ")))
  extracol <- setdiff(names(siteInfo), expectedcols)
  if(length(extracol) > 0) stop(paste("unexpected columns in siteInfo file:", paste0("'", extracol, "'", collapse=", ")))
  
  # convert dates to Date
  siteInfo <- mutate(
    siteInfo,
    date.start = as.Date(date.start, format='%Y-%m-%d'),
    date.end = as.Date(date.end, format='%Y-%m-%d')
  )
  
  # identify the constituent and flow variable named in siteInfo and inputs
  obsVars <- unique(siteInfo$constituent)
  flow <- inputs$discharge
  constitsSiteinfo <- setdiff(obsVars, flow)
  constitsInputs <- inputs$constituents
  constitsFinal <- intersect(constitsSiteinfo, constitsInputs)
  
  # flag the flow rows
  siteInfo <- mutate(siteInfo, is.flow = constituent == flow)
  
  # check inputs for errors
  if(!flow %in% obsVars) {
    stop('siteInfo constituents must include discharge as named in inputs yaml (', inputs$discharge, ')')
  }
  if(length(inpOnly <- setdiff(constitsInputs, constitsSiteinfo)) > 0) {
    warning('some constituents named in inputs yaml are not in siteInfo: ', paste0(inpOnly, collapse=', '))
  }
  
  # subset siteInfo if specified in inputs
  if(length(siOnly <- setdiff(constitsSiteinfo, constitsInputs)) > 0) {
    message('subsetting siteInfo to only those constituents named in inputs yaml\n',
            '  removing: ', paste0(siOnly, collapse=', ', '\n'),
            '  keeping: ', paste0(constitsFinal, collapse=', '))
    siteInfo <- filter(siteInfo, is.flow | constituent %in% constitsFinal)
  }
  
  # require that all specified constituent and discharge folders exist
  allFolders <- file.path(inputs$inputFolder, c(constitsFinal, flow))
  if(!all(dir.exists(allFolders))) {
    stop("Input or constituent folder does not exist")
  }
  
  # compute expected filenames based on siteInfo
  siteInfo <- siteInfo %>%
    mutate(filepath = file.path(inputs$inputFolder, constituent, paste0(site.id, '.csv')))
  
  # check that all named files exist; remove those that don't exist (#142)
  if(length(missingFiles <- siteInfo$filepath[!file.exists(siteInfo$filepath)]) > 0) {
    warning("omitting these missing files from the analysis:\n", paste0('  ', missingFiles, collapse='\n'))
    siteInfo <- filter(siteInfo, !(filepath %in% missingFiles))
  }
  
  return(siteInfo)
} 

#' Create a data.frame of site and file info where each row refers to a single 
#' site-constituent-discharge combo to be modeled
#' 
#' @param siteInfo data.frame of site information, augmented by combineSpecs()
matchFiles <- function(siteInfo) {
  
  # get the unique constituents in siteInfo
  siteConstPairs <- siteInfo %>%
    filter(!is.flow) %>%
    select(site.id, constituent) %>%
    distinct
  
  siteFileSets <- bind_rows(lapply(seq_len(nrow(siteConstPairs)), function(siteConstRow) {
    site.id <- siteConstPairs[siteConstRow, 'site.id']
    constituent <- siteConstPairs[siteConstRow, 'constituent']
    
    # find the row with constituent site/file info
    constFile <- siteConstPairs[siteConstRow,] %>%
      left_join(siteInfo, by=c('site.id','constituent'))
    if(nrow(constFile) != 1) {
      stop('need 1 site-constituent file, found:\n', paste0('  ', constFile$filepath, collapse='\n'))
    }
    
    # find the row with flow site/file info
    flowFile <- select(constFile, matching.site) %>%
      left_join(filter(siteInfo, is.flow), by='matching.site')
    if(nrow(flowFile) != 1) {
      stop('need 1 flow file for ', site.id, ' ', constituent, ', found:\n', paste0('  ', flowFile$filepath, collapse='\n'))
    }
    
    # return their joined 1-row df
    full_join(constFile, flowFile, suffix=c('.CONC', '.FLOW'), by='matching.site') %>%
      select(matching.site, everything())
  }))
  
  return(siteFileSets)
}

#' Summarize the input data for all models
#' 
#' Compute and save info on the site, constituent, and input datasets.
#' Compute num.censored specially here because we're using rloadest format
#' for censored data (smwrQW format) and will eventually have something
#' simpler in place for loadflex, at which point we'll add that to
#' summarizeInputs.
#' @param siteMeta metadata for this site and constituent
#' @param siteConstit data.frame of consitutent concentration & discharge observations
#' @param siteQ data.frame of dates and discharges
#' @param loadflexVersion version of loadflex being used
#' @param batchStartTime datetime this run was started
summarizeBatchInputs <- function(siteMeta, siteConstit, siteQ, loadflexVersion, batchStartTime) {
  inputMetrics <- summarizeInputs(siteMeta, fitdat=siteConstit, estdat=siteQ)
  inputMetrics$fitdat.num.censored <- length(which(siteConstit[['status']] == 2))
  inputMetrics$estdat.num.censored <- NULL # assuming no censoring in Q
  inputMetrics <- inputMetrics %>%
    mutate(
      loadflex.version = loadflexVersion, 
      run.date = batchStartTime)
}

#' Produce metrics describing the performance of each model
#' 
#' Describe model performance in a 1-row data.frame
#' 
#' @param allModels list of fitted model objects
#' @param siteMeta metadata for this site and constituent
#' @param loadflexVersion version of loadflex being used
#' @param batchStartTime datetime this run was started
summarizeMetrics <- function(allModels, siteMeta, loadflexVersion, batchStartTime) {
  summarize_model <- function(load.model) {
    switch( 
      class(load.model), # slightly different call for each model type
      loadReg2 = summarizeModel(load.model),
      loadComp = summarizeModel(load.model, newdata=siteQ, irregular.timesteps.ok=TRUE),
      loadInterp = summarizeModel(load.model, irregular.timesteps.ok=TRUE),
      loadBeale = data_frame(site.id=siteMeta@site.id, constituent=siteMeta@constituent, nstrata=load.model$nstrata)
    )
  }
  metrics <- bind_cols(
    as_data_frame(summarize_model(allModels[[1]])[1:2]), # site.id and constituent columns just once
    lapply(names(allModels), function(mod) { # model-specific columns
      modSum <- summarize_model(allModels[[mod]])
      modSum[-(1:2)] %>% # don't duplicate the site.id and constituent columns
        setNames(paste0(mod, '.', names(.))) # add the "REG.", "CMP.", etc. prefix
    })) %>%
    mutate(
      loadflex.version = loadflexVersion, 
      run.date = batchStartTime)
  return(metrics)
}

#' Prepare estimation data for a loadReg2 regression model
#'
#' Compute inputs for fully specified model equations (not model(7) but lnQ,
#' DECTIME, etc.), possibly with fixed DECTIME
#'
#' @param mod a fitted loadReg2 model
#' @param estdat a data.frame of inputs for load prediction
#' @param regBaseYear an integer water year to which predictions should be fixed
#'   (e.g., 2006 to fix all dates to 2006-04-01) or NA to leave dates unfixed
prepareRegEstdat <- function(mod, estdat, regBaseYear) {
  Qadj <- getFittedModel(mod)$Qadj
  Tadj <- getFittedModel(mod)$Tadj
  RLestdat <- rloadest:::setXLDat(data=estdat, flow=qColName, dates=dateColName, Qadj=Qadj, Tadj=Tadj, model.no=9) %>%
    as.data.frame() %>%
    bind_cols(estdat)
  if(!is.na(regBaseYear)) {
    estdatF <- estdat[1,]
    estdatF$fixed.date <- as.Date(sprintf("%s-04-01", regBaseYear))
    RLestdatF <- rloadest:::setXLDat(data=estdatF, flow=qColName, dates="fixed.date", Qadj=Qadj, Tadj=Tadj, model.no=9)
    RLestdat$DECTIME <- RLestdatF[[1,'DECTIME']]
    RLestdat$DECTIME2 <- RLestdatF[[1,'DECTIME2']]
    # leave sin.DECTIME and cos.DECTIME alone because we want to keep seasonality
  }
  return(RLestdat)
}

#' Produce daily load estimates
#' 
#' @param allModels list of fitted model objects
#' @param siteQ data.frame of dates and discharges
#' @param conv.load.rate multiplier for converting loads as predicted from models to loads requested by batch user
#' @param regBaseYear an integer water year to which predictions should be fixed
#'   (e.g., 2006 to fix all dates to 2006-04-01) or NA to leave dates unfixed
summarizeDaily <- function(allModels, siteQ, conv.load.rate, regBaseYear) {
  predsLoad <- lapply(allModels, function(mod) {
    (if(is(mod, 'loadReg2')) {
      RLsiteQ <- prepareRegEstdat(mod, siteQ, regBaseYear)
      predictSolute(mod, "flux", RLsiteQ, se.pred=TRUE, date=TRUE)
    } else if(is(mod, 'loadComp')) {
      suppressWarnings(predictSolute(mod, "flux", siteQ, se.pred=FALSE, date=TRUE)) %>%
        mutate(se.pred=NA) # avoiding se.pred not because we can't compute it but because it's slow to compute the MSE and we won't be using se.pred
    } else if(is(mod, 'loadBeale')) {
      data_frame(date=siteQ[[dateColName]], fit=NA, se.pred=NA)
    } else {
      predictSolute(mod, "flux", siteQ, se.pred=TRUE, date=TRUE)
    }) %>%
      mutate(
        fit = fit * conv.load.rate,
        se.pred = se.pred * conv.load.rate
      )
  })
  return(predsLoad)
}

#' Produce monthly load estimates
#' 
#' @param allModels list of fitted model objects
#' @param predsLoad list of data.frames of daily predictions
#' @param inputs list of configuration inputs
#' @param siteQ data.frame of dates and discharges
#' @param conv.load.rate multiplier for converting loads as predicted from models to loads requested by batch user
#' @param loadflexVersion version of loadflex being used
#' @param batchStartTime datetime this run was started
#' @param regBaseYear an integer water year to which predictions should be fixed
#'   (e.g., 2006 to fix all dates to 2006-04-01) or NA to leave dates unfixed
summarizeMonthly <- function(allModels, predsLoad, inputs, siteQ, conv.load.rate, regBaseYear, loadflexVersion, batchStartTime) {
  message(" * generating monthly mean load estimates...")
  monthlySummary <- bind_rows(lapply(names(predsLoad), function(mod) {
    message(paste0('   ', mod, '...'), appendLF = FALSE)
    if(is(allModels[[mod]], 'loadReg2')) {
      siteQ %>%
        mutate(Month=format(siteQ[[dateColName]], '%Y-%m')) %>%
        group_by(Month) %>%
        do({
          # months with a lot of NaNs in se.pred really slow rloadest down. 
          # skip months exceeding the criterion (inputs$regMaxNaNsPerMonth)
          siteMonthQ <- .
          dailyPreds <- predsLoad[[mod]] %>% 
            mutate(Month=format(predsLoad[[mod]][[dateColName]], '%Y-%m')) %>%
            filter(Month == siteMonthQ$Month[1])
          if(length(which(!is.finite(dailyPreds$se.pred))) > inputs$regMaxNaNsPerMonth) {
            message('skipping NaN-riddled ', as.character(siteMonthQ$Month[1]), '...', appendLF = FALSE)
            data_frame(Flux=NaN, SEP=NaN, Ndays=as.numeric(NA))
          } else {
            # filter to non-NA se.preds
            siteMonthQ <- filter(siteMonthQ, is.finite(dailyPreds$se.pred)) %>%
              select(-Month)
            # compute inputs for fully specified model equations (not model(7)
            # but lnQ, DECTIME, etc.), possibly with fixed DECTIME
            RLsiteMonthQ <- prepareRegEstdat(allModels[[mod]], siteMonthQ, regBaseYear)
            # pass filtered data for fully-specified equation to model
            predLoad(getFittedModel(allModels[[mod]]), newdata=RLsiteMonthQ, by='total', allow.incomplete=TRUE)
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
        select(Month, Flux_Rate, SE, n, CI_lower, CI_upper, model)
    } else if(is(allModels[[mod]], 'loadBeale')) {
      data_frame(
        Month=unique(format(siteQ[[dateColName]], '%Y-%m')),
        Flux_Rate = NA,
        SE = NA,
        n = NA,
        CI_lower = NA,
        CI_upper = NA,
        model=mod)
    } else {
      suppressWarnings(loadflex:::aggregateSolute(
        predsLoad[[mod]], siteMeta, agg.by="month", format='flux rate')) %>%
        mutate(
          Month = as.character(Month),
          SE = NA,
          CI_lower = NA,
          CI_upper = NA,
          model=mod)
    }
  })) %>%
    mutate(site.id=getInfo(siteMeta, 'site.id'), constituent=getInfo(siteMeta, 'constituent')) %>%
    select(site.id, constituent, model, everything())
  col_order <- c(
    'site.id', 'constituent', 'Month',
    sapply(names(allModels), function(mod) paste0(mod, '.', c('n', 'Flux_Rate', 'SE', 'CI_lower', 'CI_upper'))))
  monthlySummary <- monthlySummary %>%
    tidyr::gather(var, val, Flux_Rate, SE, n, CI_lower, CI_upper) %>%
    mutate(var=ordered(paste0(model, '.', var))) %>%
    select(-model) %>%
    tidyr::spread(var, val) %>%
    select_(.dots=col_order)
  monthlySummary <- monthlySummary %>%
    mutate(
      loadflex.version = loadflexVersion, 
      run.date = batchStartTime)
  message('done!')
  return(monthlySummary)
}

#' Produce seasonal load estimates
#' 
#' @param allModels list of fitted model objects
#' @param predsLoad list of data.frames of daily predictions
#' @param inputs list of configuration inputs
#' @param siteQ data.frame of dates and discharges
#' @param conv.load.rate multiplier for converting loads as predicted from models to loads requested by batch user
#' @param loadflexVersion version of loadflex being used
#' @param batchStartTime datetime this run was started
summarizeSeasonal <- function(allModels, predsLoad, inputs, siteQ, conv.load.rate, regBaseYear, loadflexVersion, batchStartTime) {
  message(" * generating seasonal mean load estimates...")
  as_season <- function(dates) {
    ends <- c('03/30','06/30','09/30','12/31')
    starts <- c('01-01','04-01','07-01','10-01')
    breaknames <- paste0(starts, '-', ends)
    seasons <- smwrBase::seasons(dates, breaks=ends, Names=starts)
    paste0(format(dates, "%Y"), '-', seasons)
  }
  seasonalSummary <- bind_rows(lapply(names(predsLoad), function(mod) {
    message(paste0('   ', mod, '...'), appendLF = FALSE)
    if(is(allModels[[mod]], 'loadReg2')) {
      siteQ %>%
        mutate(Season=as_season(siteQ[[dateColName]])) %>%
        group_by(Season) %>%
        do({
          # months with a lot of NaNs in se.pred really slow rloadest down. 
          # skip months exceeding the criterion (inputs$regMaxNaNsPerSeason)
          siteSeasonQ <- .
          dailyPreds <- predsLoad[[mod]] %>% 
            mutate(Season=as_season(predsLoad[[mod]][[dateColName]])) %>%
            filter(Season == siteSeasonQ$Season[1])
          if(length(which(!is.finite(dailyPreds$se.pred))) > inputs$regMaxNaNsPerSeason) {
            message('skipping NaN-riddled ', as.character(siteSeasonQ$Season[1]), '...', appendLF = FALSE)
            data_frame(Flux=NaN, SEP=NaN, Ndays=as.numeric(NA))
          } else {
            siteSeasonQ <- filter(siteSeasonQ, is.finite(dailyPreds$se.pred)) %>%
              select(-Season)
            # compute inputs for fully specified model equations (not model(7)
            # but lnQ, DECTIME, etc.), possibly with fixed DECTIME
            RLsiteSeasonQ <- prepareRegEstdat(allModels[[mod]], siteSeasonQ, regBaseYear)
            predLoad(getFittedModel(allModels[[mod]]), newdata=RLsiteSeasonQ, by='total', allow.incomplete=TRUE)
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
        select(Season, Flux_Rate, SE, n, CI_lower, CI_upper, model)
    } else if(is(allModels[[mod]], 'loadBeale')) {
      data_frame(
        Season=unique(as_season(siteQ[[dateColName]])),
        Flux_Rate = NA,
        SE = NA,
        n = NA,
        CI_lower = NA,
        CI_upper = NA,
        model=mod)
    } else {
      predsSeason <- predsLoad[[mod]] %>% mutate(Season=as_season(predsLoad[[mod]][[dateColName]]))
      suppressWarnings(loadflex:::aggregateSolute(
        predsSeason,
        siteMeta, custom=predsSeason['Season'], agg.by="Season", format='flux rate')) %>%
        mutate(
          Season = as.character(Season),
          SE = NA,
          CI_lower = NA,
          CI_upper = NA,
          model=mod) %>%
        as_data_frame()
    }
  })) %>%
    mutate(site.id=getInfo(siteMeta, 'site.id'), constituent=getInfo(siteMeta, 'constituent')) %>%
    select(site.id, constituent, model, everything())
  col_order <- c(
    'site.id', 'constituent', 'Season',
    sapply(names(allModels), function(mod) paste0(mod, '.', c('n', 'Flux_Rate', 'SE', 'CI_lower', 'CI_upper'))))
  seasonalSummary <- seasonalSummary %>%
    tidyr::gather(var, val, Flux_Rate, SE, n, CI_lower, CI_upper) %>%
    mutate(var=ordered(paste0(model, '.', var))) %>%
    select(-model) %>%
    tidyr::spread(var, val) %>%
    select_(.dots=col_order)
  seasonalSummary <- seasonalSummary %>%
    mutate(
      loadflex.version = loadflexVersion, 
      run.date = batchStartTime)
  message('done!')
  return(seasonalSummary)
}
#' Produce annual load estimates
#' 
#' @param allModels list of fitted model objects
#' @param predsLoad list of data.frames of daily predictions
#' @param inputs list of configuration inputs
#' @param siteQ data.frame of dates and discharges
#' @param conv.load.rate multiplier for converting loads as predicted from models to loads requested by batch user
#' @param loadflexVersion version of loadflex being used
#' @param batchStartTime datetime this run was started
summarizeAnnual <- function(allModels, predsLoad, inputs, siteQ, conv.load.rate, regBaseYear, loadflexVersion, batchStartTime) {
  message(" * generating annual mean load estimates...")
  annualSummary <- bind_rows(lapply(names(predsLoad), function(mod) {
    message(paste0('   ', mod, '...'), appendLF = FALSE)
    if(is(allModels[[mod]], 'loadReg2')) {
      siteQ %>%
        mutate(Water_Year=as.character(smwrBase::waterYear(siteQ[[dateColName]]))) %>%
        group_by(Water_Year) %>%
        do({
          # years with a lot of NaNs in se.pred really slow rloadest down. 
          # skip years exceeding the criterion (inputs$regMaxNaNsPerYear)
          siteYearQ <- .
          dailyPreds <- predsLoad[[mod]] %>% 
            mutate(Water_Year=smwrBase::waterYear(predsLoad[[mod]][[dateColName]])) %>%
            filter(Water_Year == siteYearQ$Water_Year[1])
          if(length(which(!is.finite(dailyPreds$se.pred))) > inputs$regMaxNaNsPerYear) {
            message('skipping NaN-riddled ', as.character(siteYearQ$Water_Year[1]), '...', appendLF = FALSE)
            data_frame(Flux=NaN, SEP=NaN, Ndays=as.numeric(NA)) # Ndays=NA ensures this year is skipped for multi-year estimate, too
          } else {
            siteYearQ <- filter(siteYearQ, is.finite(dailyPreds$se.pred)) %>%
              select(-Water_Year)
            # compute inputs for fully specified model equations (not model(7)
            # but lnQ, DECTIME, etc.), possibly with fixed DECTIME
            RLsiteYearQ <- prepareRegEstdat(allModels[[mod]], siteYearQ, regBaseYear)
            predLoad(getFittedModel(allModels[[mod]]), newdata=RLsiteYearQ, by='water year', allow.incomplete=TRUE)
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
    } else if(is(allModels[[mod]], 'loadBeale')) {
      data_frame(
        Water_Year=as.character(unique(smwrBase::waterYear(siteQ[[dateColName]]))),
        Flux_Rate = NA,
        SE = NA,
        n = NA,
        CI_lower = NA,
        CI_upper = NA,
        model=mod)
    } else {
      suppressWarnings(loadflex:::aggregateSolute(
        predsLoad[[mod]], siteMeta, agg.by="water year", format='flux rate')) %>%
        mutate(
          Water_Year = as.character(Water_Year),
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
  return(annualSummary)
}

#' Produce annual load estimates
#' 
#' @param allModels list of fitted model objects
#' @param predsLoad list of data.frames of daily predictions
#' @param inputs list of configuration inputs
#' @param siteQ data.frame of dates and discharges
#' @param conv.load.rate multiplier for converting loads as predicted from models to loads requested by batch user
#' @param loadflexVersion version of loadflex being used
#' @param batchStartTime datetime this run was started
summarizeMultiYear <- function(allModels, predsLoad, annualSummary, inputs, siteQ, conv.load.rate, regBaseYear, loadflexVersion, batchStartTime) {
  message(" * generating multi-year mean load estimates...", appendLF = FALSE)
  multiYearSummary <- bind_rows(lapply(names(predsLoad), function(mod) {
    message(paste0(mod, '...'), appendLF=FALSE)
    completeWaterYears <- annualSummary %>%
      filter(
        !is.na(.[[paste0(mod,'.n')]]), # .n=NA if water year was skipped for inputs$regMaxNaNsPerYear
        .[[paste0(mod,'.n')]] >= inputs$minDaysPerYear) %>%
      .$Water_Year
    completeSiteQ <- mutate(siteQ, Water_Year = smwrBase::waterYear(date)) %>%
      filter(Water_Year %in% completeWaterYears)
    (if(is(allModels[[mod]], 'loadReg2')) {
      # compute inputs for fully specified model equations (not model(7)
      # but lnQ, DECTIME, etc.), possibly with fixed DECTIME
      RLcompleteSiteQ <- prepareRegEstdat(allModels[[mod]], completeSiteQ, regBaseYear)
      predLoad(getFittedModel(allModels[[mod]]), newdata=RLcompleteSiteQ, by='total', allow.incomplete=TRUE) %>%
        mutate(
          Flux_Rate = Flux * conv.load.rate,
          SE = SEP * conv.load.rate,
          CI_lower = L95 * conv.load.rate,
          CI_upper = U95 * conv.load.rate,
          years.record = length(unique(annualSummary$Water_Year)),            
          years.complete = length(unique(completeWaterYears))
        ) %>%
        select(Flux_Rate, SE, CI_lower, CI_upper, years.record, years.complete)
    } else if(is(allModels[[mod]], 'loadBeale')) {
      data_frame(
        Flux_Rate = allModels[[mod]]$rload_kg_y,
        SE = allModels[[mod]]$serload_kg_y
      ) %>%
        mutate(
          CI_lower = Flux_Rate - 1.96*SE,
          CI_upper = Flux_Rate + 1.96*SE,
          years.record = length(unique(annualSummary$Water_Year)),            
          years.complete = length(unique(completeWaterYears)))
    } else {
      suppressWarnings(loadflex:::aggregateSolute(
        predsLoad[[mod]], siteMeta, agg.by="mean water year", 
        format='flux rate', min.n=inputs$minDaysPerYear, ci.agg=FALSE, se.agg=FALSE))
    }) %>%
      mutate(model=mod)
  })) %>%
    mutate(site.id=getInfo(siteMeta, 'site.id'), constituent=getInfo(siteMeta, 'constituent')) %>%
    select(site.id, constituent, model, everything())
  multiYearSummary <- reshape(
    multiYearSummary, idvar = c('site.id','constituent'), direction = "wide", 
    v.names = c("Flux_Rate", "SE", "CI_lower", "CI_upper",'years.record','years.complete'), timevar = "model")
  multiYearSummary <- setNames( # replace the ".REG" or ".CMP" suffixes with "REG.", "CMP.", prefixes
    multiYearSummary, 
    sub(pattern=sprintf('(.*)\\.(%s)', paste0(names(allModels), collapse='|')), 
        replacement='\\2.\\1', names(multiYearSummary)))
  multiYearSummary <- multiYearSummary %>%
    mutate(
      loadflex.version = loadflexVersion, 
      run.date = batchStartTime)
  message('done!')
  return(multiYearSummary)
}

#' Combine 1-row summary files into a single multi-row file
#' 
#' @param csvType the name of the summary type to combine
#' @param constitSiteInfo the metadata table linking constituent and discharge files
#' @param outputFolder the folder where output should be written
summarizeCsvs <- function(csvType=c('inputs','annual','multiYear', 'modelMetrics'), 
                          constitSiteInfo, outputFolder) {
  csvType <- match.arg(csvType)
  
  # read and combine the 1-row data files for all sites and constituents
  allCsvs <- bind_rows(lapply(seq_len(nrow(constitSiteInfo)), function(i) {
    matchingSite <- constitSiteInfo$matching.site[i] # this is how we'll name the output files
    constitName <- constitSiteInfo$constituent.CONC[i]
    csvFile <- file.path(outputFolder, constitName, csvType, paste0(matchingSite, '.csv'))
    tryCatch({
      read.csv(csvFile, header=TRUE, stringsAsFactors=FALSE)
    }, error=function(e) {
      message(paste0('could not read ', csvFile), call.=FALSE)
      NULL
    })
  }))
  
  # write separate csvs for each constituent
  constits <- unique(allCsvs$constituent)
  sapply(constits, function(con) {
    writeFile = file.path(outputFolder, con, paste0(con,"_", csvType,".csv"))
    write.csv(filter(allCsvs, constituent == con), 
              file = writeFile, 
              row.names = FALSE)
    message('the ', con, ' ', csvType, ' summary has been written to ', writeFile)
  })
  
  return(allCsvs)
}

#' Combine 1-row summary files into a single multi-row file
#' 
#' @param constitSiteInfo the metadata table linking constituent and discharge files
#' @param outputFolder the folder where output should be written
summarizePlots <- function(constitSiteInfo, outputFolder) {
  constits <- unique(constitSiteInfo$constituent.CONC)
  for(constitName in constits) {
    constitFileSets <- dplyr::filter(constitSiteInfo, constituent.CONC == constitName)
    # read and combine the 1-site, 1 constit pdf files for all sites with this constituent
    allPlotFiles <- unlist(lapply(seq_len(nrow(constitFileSets)), function(i) {
      matchingSite <- constitSiteInfo$matching.site[i] # this is how we named the output files
      plotFile <- file.path(outputFolder, constitName, 'plots', paste0(matchingSite, '.pdf'))
      if(file.exists(plotFile)) {
        plotFile
      } else {
        message(paste0('could not read ', plotFile))
        NULL
      }
    }))
    if(length(allPlotFiles) > 0) {
      if(require(plotflow)) {
        ## paste the paths to pdfs together in one string w/ spaces
        plotflow:::mergePDF(
          in.file=paste(allPlotFiles, collapse=" "),
          file=paste0(constitName, "_plots.pdf")
        )
      } else {
        message("to combine pdfs, install the plotflow package with devtools::install_github('trinker/plotflow')")
      }
    } else {
      message("no ", constitName, " pdfs to combine")
    }
  }
}

#' Create plots for all models for a single site/constituent pair
#' 
#' @param loadModels a list of load models
#' @param estdata data.frame of estimation data (dates and discharges)
#' @param siteMeta loadflex metadata object
writePDFreport <- function(loadModels, fitdat, censdat, estdat, siteMeta, loadflexVersion, batchStartTime) {
  
  # make plots. the first page is redundant across models
  modelNames <- data_frame(
    short = c("RL5", "RL7", "INT", "CMP"),
    long = c("Regression Model L5 (rloadest 5 parameter)",
             "Regression Model L7 (rloadest 7 parameter)",
             "Interpolation Model (rectangular)",
             "Composite Model (rloadest + interpolation)"))
  
  # page 1: input data with censoring
  eList <- suppressWarnings( # ANA example data: This program requires at least 30 data points. Rolling means will not be calculated.
    convertToEGRET(meta = siteMeta, data=censdat, newdata = estdat))
  plotEGRET("multiPlotDataOverview", eList=eList)
  title(main="Input Data", line=-1, adj=0, outer=TRUE)
  title(main=sprintf("%s-%s", siteMeta@site.id, 'ALL'), line=-1, adj=1, outer=TRUE)
  mtext(text=sprintf('loadflex version %s', loadflexVersion), side=1, line=-1, adj=0, outer=TRUE, font=3)
  mtext(text=sprintf('run at %s', batchStartTime), side=1, line=-1, adj=1, outer=TRUE, font=3)
  
  # create expanded siteQ data for rloadest models
  rlmodids <- which(sapply(allModels, is, 'loadReg2'))
  if(length(rlmodids) > 0) {
    loadModel <- allModels[[rlmodids[1]]]
    Qadj <- getFittedModel(loadModel)$Qadj
    Tadj <- getFittedModel(loadModel)$Tadj
    RLestdat <- rloadest:::setXLDat(data=estdat, flow=qColName, dates=dateColName, Qadj=Qadj, Tadj=Tadj, model.no=9) %>%
      as.data.frame() %>%
      bind_cols(estdat)
  }

  for(m in which(!sapply(allModels, is, 'loadBeale'))) {
    loadModel <- loadModels[[m]]
    loadModel@metadata <- siteMeta
    modelShort <- modelNames$short[modelNames$short == names(loadModels)[m]]
    modelLong <- modelNames$long[modelNames$short == names(loadModels)[m]]
    if(is(loadModel, 'loadReg2')) {
      eList <- suppressWarnings( # ANA example data: This program requires at least 30 data points. Rolling means will not be calculated.
        convertToEGRET(load.model=loadModel, newdata=RLestdat))
    } else {
      eList <- suppressWarnings( # CMP: Uncertainty estimates are unavailable. Proceeding with NAs
        convertToEGRET(load.model=loadModel, newdata=estdat))
    }
    
    # pages 2,4,6,8
    par(mfrow=c(2,1), oma=c(0,0,1,0))
    plotEGRET("plotConcTimeDaily", eList=eList, mgp = c(4,1,0))
    title(main=paste0("Predictions"), line=0, adj=0, outer=TRUE)
    title(main=sprintf("%s-%s-1", siteMeta@site.id, modelShort), line=0, adj=1, outer=TRUE)
    plotEGRET("plotFluxTimeDaily", eList=eList, mgp = c(4,1,0))
    
    # pages 3,5,7,9
    if(is(loadModel, 'loadReg2')) {
      plotEGRET("fluxBiasMulti", eList=eList, moreTitle = paste0(modelLong, '; '))
      headerLine <- -1
    } else {
      par(oma=c(0,0,1,0))
      layout(matrix(c(0,0,0,0, 0,1,2,0, 0,3,4,0, 0,0,0,0), nrow=4, ncol=4, byrow=TRUE), widths=c(0.15,0.3,0.3,0.15))
      plotEGRET("boxConcThree", eList=eList, printTitle = FALSE, tinyPlot = TRUE)
      plotEGRET("plotConcPred", eList=eList, printTitle = FALSE, tinyPlot = TRUE, randomCensored = FALSE)
      plotEGRET("boxQTwice", eList=eList, printTitle = FALSE, tinyPlot = TRUE)
      plotEGRET("plotFluxPred", eList=eList, printTitle = FALSE, tinyPlot = TRUE, randomCensored = FALSE)
      fluxBias <- as.numeric(EGRET::fluxBiasStat(eList$Sample)[3])
      title <- paste(
        eList$INFO$shortName, ", ", eList$INFO$paramShortName, 
        "\nModel is ", modelLong, "; Flux Bias Statistic = ", format(fluxBias, digits=3), 
        sep = "")
      par(cex=1, cex.main=1.2)
      mtext(title, outer = TRUE, line = -18, font = 1.8)
      headerLine <- -0.5
    }
    title(main=paste0("Diagnostics"), line=headerLine, adj=0, outer=TRUE, cex=1)
    title(main=sprintf("%s-%s-2", siteMeta@site.id, modelShort), line=headerLine, adj=1, outer=TRUE, cex=1)
  }
}
