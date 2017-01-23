#TODO: tests!

library(dplyr)
library(loadflex)
library(tools)
library(rloadest)

#run loadflex over multiple sites

#------------------User Inputs--------------------#

#TODO: implement this
outputFormat <- "simple" #or "complex"

#input constituents
#script will look for folders with this name, and use in site metadata
constituents <- c("NO3", "PT")
loadUnits <- "kg"
loadRateUnits <- "kg/d"

inputFolder <- "three_ANA_sites" #folder containing all input subfolders
outputFolder <- "./output"
dischargeFolder <- "Q" #subfolder of inputFolder containing discharge measurements for predictions
siteInfo <- "siteInfo.csv" #also inside inputFolder, data frame of site info

#-------------------------Check files, set up directories-----------------------# 
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

fileDF <- makeFileDF(inputFolder, constits = constituents, discharge.folder = dischargeFolder)
allSiteInfo <- read.csv(file.path(inputFolder, siteInfo), stringsAsFactors = FALSE)

#setup output directories
nConstits <- length(constituents)
outConstit <- file.path(rep(outputFolder, nConstits), constituents)
sapply(outConstit, dir.create, recursive = TRUE, showWarnings = FALSE)
outTemporal <- file.path(rep(outConstit,3), c(rep("inputs", nConstits), rep("annual", nConstits), rep("multiYear", nConstits)))
sapply(outTemporal, dir.create, showWarnings = FALSE)

#-----------------loadflex--------------#

siteSummaries <- data.frame()
modelSummaries <- data.frame()
allModels <- list()
annuals <- data.frame()

#loop over unique sites
for(i in 1:nrow(fileDF)) {
  message(paste('processing constituent file', fileDF$constitFile[i], '\n'))
  
  #TODO: use siteInfo to get column names, fail if they aren't the same as directories
  
  #read in appropriate files
  siteQ <- read.csv(fileDF$qFile[i], stringsAsFactors = FALSE)
  siteConstit <- read.csv(fileDF$constitFile[i], stringsAsFactors = FALSE)
  
  #needs to convert dates 
  siteConstit$date <- as.Date(siteConstit$date)
  siteQ$date <- as.Date(siteQ$date)
  
  #pull out appropriate rows of allSiteInfo for Q and constit
  #need to extract constit and stations from file paths, so we know 
  #what row of site info to look at
  #TODO: deal with different discharge/consituent drainage areas here?
  constitStation <- basename(file_path_sans_ext(fileDF$constitFile[i])) 
  constitName <- basename(dirname(fileDF$constitFile[i]))
  
  constitSiteInfo <- filter(allSiteInfo, file == constitStation, constituent == constitName)
  qSiteInfo <- filter(allSiteInfo, file == constitStation, constituent == 'Q')
  
  #create metadata
  #not sure units etc are following the correct format
  #currently depending on consistent column positions to get names
  constitColName <- names(siteConstit)[3]
  qColName <- names(siteConstit)[2]
  dateColName <- names(siteConstit)[1]
  siteMeta <- metadata(constituent = constitColName, flow = qColName, dates = dateColName,
                       conc.units = constitSiteInfo$constit_units, flow.units = qSiteInfo$constit_units, load.units = loadUnits,
                       load.rate.units = loadRateUnits, station = constitStation, 
                       consti.name = "test")
    
  #TODO: site metrics
  siteMetrics <- summarizeSite(constitSiteInfo, siteConstit)
  # compute and save info on the site, constituent, and input datasets (we'll
  # recombine in the next loop)
  inputMetrics <- summarizeInputs(siteMeta, fitdat=siteConstit, estdat=siteQ)
  write.csv(inputMetrics, file.path(outputFolder, constitName, "inputs", paste0(constitStation, '.csv')), row.names=FALSE)
  
  #fit models
  #TODO: decide on standard column names?  user input timestep above?
  loadRegFormula <- formula(paste(constitColName,"~model(7)"))
  rloadest5param <- loadReg2(loadReg(loadRegFormula, data = siteConstit[1:3], 
                                     flow = qColName, dates = dateColName, time.step = "day",
                                     flow.units = getInfo(siteMeta, 'flow.units', unit.format = "rloadest"), 
                                     conc.units = getInfo(siteMeta, 'conc.units', unit.format = "rloadest"),
                                     load.units = getInfo(siteMeta, 'load.units')))
  
  interpRect <- loadInterp(interp.format = "conc", interp.function = rectangularInterpolation,
                           data = siteConstit, metadata = siteMeta)
  comp <- loadComp(reg.model = rloadest5param, interp.format = "conc", interp.function = rectangularInterpolation, 
                   interp.data = siteConstit)
  
  #list of all model objects
  allModels[[constitStation]] <- list(comp = comp, interpRect = interpRect, 
                                      rloadest5param = rloadest5param)
  
  #make predictions
  pred_rload <- predictSolute(rloadest5param, "flux", siteQ, 
                              se.pred = TRUE, date = TRUE)
  pred_interp <- predictSolute(interpRect, "flux", siteQ, 
                               se.pred = TRUE, date = TRUE)
  pred_comp <- predictSolute(comp, "flux", siteQ, se.pred = TRUE,
                             date = TRUE)
  
  #TODO: model metrics
  annualPreds <- bind_rows(
    summarizePreds(pred_rload, siteMeta, "total", model.name = "rloadest"),
    summarizePreds(pred_interp, siteMeta, "total", model.name = "interpolation"),
    summarizePreds(pred_comp, siteMeta, "total", model.name = "composite"))
  
  #TODO: plots
  
  
  #TODO: recombine into single dfs
   siteSummaries <- bind_rows(siteSummaries, siteMetrics)
   annuals <- bind_rows(annuals, annualSite)
   
   #TODO: verbose option to print output?
   write.csv(x = siteMetrics, file = file.path(outputFolder, constitName, 
                                               "multiYear", paste0(constitStation, '.csv')))
   write.csv(x = annuals, file = file.path(outputFolder, constitName, 
                                           "annual", paste0(constitStation, '.csv')))
   message(paste('Finished processing constituent file', fileDF$constitFile[i], '\n'))
}

