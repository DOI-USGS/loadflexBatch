#read-in function

# make a data.frame describing the data files that exist in this directory for
# the given constituents
makeFileDF <- function(input.folder, constits, discharge.folder) {
  #TODO:check that all constituent files have matching discharge records
  #df of corresponding files
  allFolders <- file.path(input.folder, c(constits, discharge.folder))
  if(!all(dir.exists(allFolders))) {
    stop("Input or constituent folder does not exist")
  }
  
  #get all constituent files
  constitFiles <- list.files(file.path(input.folder, constits), full.names = TRUE)
  constitNameOnly <- basename(constitFiles)
  qFiles <- file.path(input.folder, dischargeFolder, constitNameOnly)
  #TODO: warning if not matching dischargeFolder, will be skipped
  #should deal with if a discharge file doesn't exist?
  
  fileDF <- data.frame(constitFile = constitFiles, qFile=qFiles, stringsAsFactors = FALSE)
  return(fileDF)
}


# recombine summaries into single dfs
summarizeCsvs <- function(csvType=c('inputs','annual','multiYear'), fileDF, outputFolder) {
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
  allCsvFile <- file.path(outputFolder, paste0(csvType, '.csv'))
  message('the summary has been written to ', allCsvFile)
  write.csv(allCsvs, allCsvFile, row.names=FALSE)
  return(allCsvs)
}

#write plots to pdfs for a single site/constituent pair
#handles "tall" preds data frame of multiple models
writePDFreport <- function(file, intdat, estdat, allPreds, meta, inputCSV, annualCSV) {
  #pdf(file, height = 11, width = 8.5)
  
  #write csv data to pretty table
  #TODO: add titles for different models
  
  # input <- tableGrob(inputCSV)
  # annual <- tableGrob(annualCSV)
  # grid.arrange(input, annual, ncol= 1)
  
  #plots
  #are some of these going to be redundant with multiple models?
  modelNames <- data.frame(short = c("rloadest", "interp", "composite"),
                           long = c("rloadest 5 parameter model",
                                    "Interpolation Model",
                                    "Composite rloadest and interpolation model"),
                           stringsAsFactors = FALSE)
  for(m in unique(allPreds$model)) {
    preds <- filter(allPreds, model == m)
    modelLong <- modelNames$long[modelNames$short == m] 
    
    #page 1
    par(omi = c(2,2,2,2))
    plotEGRET("multiPlotDataOverview", intdat, estdat, preds, meta)
    title(paste("Input data:", getInfo(siteMeta, "site.id"), modelLong))
    #page 2
    par(mfrow=c(2,1))
    plotEGRET("plotConcTimeDaily", intdat, estdat, preds, meta)
    title(paste("Predictions:", getInfo(siteMeta, "site.id"), modelLong), line = 6)
    plotEGRET("plotFluxTimeDaily", intdat, estdat, preds, meta)
    
    #page 3
    par(mfrow=c(1,1))
    plotEGRET("fluxBiasMulti", intdat, estdat, preds, meta, moreTitle = modelLong)
    title(paste("Diagnostics:", getInfo(siteMeta, "site.id"), modelLong), line = 3)
  }
  
  #dev.off()
}