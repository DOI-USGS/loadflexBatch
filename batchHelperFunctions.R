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

#' Limit the time spent on a function
#' 
#' From https://stackoverflow.com/questions/7891073/time-out-an-r-command-via-something-like-try
#' 
#' @param expr the expression to evaluate
#' @param cpu cpu argument as in setTimeLimit
#' @param elapsed argument as in setTimeLimit
tryWithTimeLimit <- function(expr, cpu = Inf, elapsed = Inf)
{
  y <- try({setTimeLimit(cpu, elapsed, transient=TRUE); expr}, silent = TRUE) 
  if(inherits(y, "try-error")) NULL else y 
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
writePDFreport <- function(loadModels, estdat, siteMeta, loadflexVersion, batchStartTime) {
  
  # make plots. the first page is redundant across models
  modelNames <- data.frame(
    short = c("RL5", "RL7", "INT", "CMP"),
    long = c("Regression Model L5 (rloadest 5 parameter)",
             "Regression Model L7 (rloadest 7 parameter)",
             "Interpolation Model (rectangular)",
             "Composite Model (rloadest + interpolation)"),
    stringsAsFactors = FALSE)
  
  # page 1: input data with censoring
  eList <- suppressWarnings( # ANA example data: This program requires at least 30 data points. Rolling means will not be calculated.
    convertToEGRET(meta = siteMeta, data=getFittingData(loadModels[[1]]), newdata = estdat))
  plotEGRET("multiPlotDataOverview", eList=eList)
  title(main="Input Data", line=-1, adj=0, outer=TRUE)
  title(main=sprintf("%s-%s", siteMeta@site.id, 'ALL'), line=-1, adj=1, outer=TRUE)
  mtext(text=sprintf('loadflex version %s', loadflexVersion), side=1, line=-1, adj=0, outer=TRUE, font=3)
  mtext(text=sprintf('run at %s', batchStartTime), side=1, line=-1, adj=1, outer=TRUE, font=3)
  
  for(m in 1:length(loadModels)) {
    loadModel <- loadModels[[m]]
    loadModel@metadata <- siteMeta
    modelShort <- modelNames$short[modelNames$short == names(loadModels)[m]]
    modelLong <- modelNames$long[modelNames$short == names(loadModels)[m]]
    eList <- suppressWarnings( # CMP: Uncertainty estimates are unavailable. Proceeding with NAs
      convertToEGRET(load.model=loadModel, newdata=estdat))
    
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
