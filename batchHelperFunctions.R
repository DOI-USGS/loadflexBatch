#' Make a data.frame describing the constituent and flow sites and the data
#' files for those variables
#' 
#' @param inputs list of user inputs as given in a yml file
combineSpecs <- function(inputs) {
  
  # read in the siteInfo file
  siteInfo <- read.csv(file.path(inputs$inputFolder, inputs$siteInfo), stringsAsFactors = FALSE)
  
  # identify the constituent and flow variable named in siteInfo and inputs
  obsVars <- unique(siteInfo$constituent)
  flow <- inputs$discharge
  constitsSiteinfo <- setdiff(obsVars, flow)
  constitsInputs <- inputs$constituents
  constitsFinal <- union(constitsSiteinfo, constitsInputs)
  
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
            '  keeping: ', paste0(si.and.inp, collapse=', '))
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

# recombine summaries into single dfs
summarizeCsvs <- function(csvType=c('inputs','annual','multiYear', 'modelMetrics'), 
                          fileDF, outputFolder) {
  csvType <- match.arg(csvType)
  allCsvs <- bind_rows(lapply(seq_len(nrow(fileDF)), function(i) {
    constitStation <- basename(file_path_sans_ext(fileDF$constitFile[i])) 
    constitName <- basename(dirname(fileDF$constitFile[i]))
    csvFile <- file.path(outputFolder, constitName, csvType, paste0(constitStation, '.csv'))
    tryCatch({
      suppressWarnings(read.csv(csvFile, header=TRUE, stringsAsFactors=FALSE))
    }, error=function(e) {
      message(paste0('could not read ', csvFile), call.=FALSE)
      NULL
    })
  }))
  
  #write separate csvs for each constituent
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

#write plots to pdfs for a single site/constituent pair
#handles "tall" preds data frame of multiple models
writePDFreport <- function(file, load.models, estdat, siteMeta) {
  #pdf(file, height = 11, width = 8.5)
  
  # make plots. the first page is redundant across models
  modelNames <- data.frame(
    short = c("rloadest", "interp", "composite"),
    long = c("rloadest 5 parameter model",
             "Interpolation Model",
             "Composite rloadest and interpolation model"),
    stringsAsFactors = FALSE)
  
  for(m in 1:length(load.models)) {
    load.model <- load.models[[m]]
    load.model@metadata <- siteMeta
    modelLong <- modelNames$long[modelNames$short == names(load.models)[m]]
    
    # page 1
    par(omi = c(2,2,2,2))
    plotEGRET("multiPlotDataOverview", load.model=load.model, newdata=estdat)
    title(paste("Input data:", getInfo(siteMeta, "site.id"), modelLong))
  
    # page 2
    par(mfrow=c(2,1))
    plotEGRET("plotConcTimeDaily", load.model=load.model, newdata=estdat, mgp = c(4,1,0))
    title(paste("Predictions:", getInfo(siteMeta, "site.id"), modelLong), line = 6)
    plotEGRET("plotFluxTimeDaily", load.model=load.model, newdata=estdat, mgp = c(4,1,0))
    
    # page 3
    par(mfrow=c(1,1))
    plotEGRET("fluxBiasMulti", load.model=load.model, newdata=estdat, moreTitle = modelLong)
    title(paste("Diagnostics:", getInfo(siteMeta, "site.id"), modelLong), line = 3)
  }
  
  #dev.off()
}
