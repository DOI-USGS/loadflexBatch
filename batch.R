# This script runs loadflex in "batch mode": It fits several load estimation 
# models for many constituents at many sites, generates and saves the 
# predictions, and produces summaries over all models, constituents, and sites.

# 1. Start by setting your working directory to the directory that contains this
# script. The input directory (three_ANA_sites) should also be within this 
# working directory.

# 2. Modify the settings within the User Inputs section to match your input data
# and preferences.

# 3. Source this script, preferably within a clean R session with no additional
# variables in the R environment.

# 4. This script will add an output directory at the path indicated by 
# outputFolder below. The current setting creates a folder called 'output' in
# the same directory as this script.

# 5. Inspect the plots and tables files in the output directory.

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

#-------------------------Load packages, check files, set up directories-----------------------# 

library(dplyr)
library(loadflex)
library(tools)
library(rloadest)
library(gridExtra)
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
outTemporal <- file.path(rep(outConstit,3), c(rep("inputs", nConstits), rep("annual", nConstits), rep("multiYear", nConstits)))
sapply(outTemporal, dir.create, showWarnings = FALSE)

#-----------------loadflex--------------#

allModels <- list() # this might get big. i'd prefer splitting & saving

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
  #need to extract constit and sites from file paths, so we know 
  #what row of site info to look at
  #TODO: deal with different discharge/consituent drainage areas here?
  constitSite <- basename(file_path_sans_ext(fileDF$constitFile[i])) 
  constitName <- basename(dirname(fileDF$constitFile[i]))
  
  constitSiteInfo <- filter(allSiteInfo, matching.site == constitSite, constituent == constitName)
  qSiteInfo <- filter(allSiteInfo, matching.site == constitSite, constituent == 'Q')
  
  #create metadata
  #not sure units etc are following the correct format
  #currently depending on consistent column positions to get names
  constitColName <- names(siteConstit)[3]
  qColName <- names(siteConstit)[2]
  dateColName <- names(siteConstit)[1]
  
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
  # recombine in the next loop)
  inputMetrics <- summarizeInputs(siteMeta, fitdat=siteConstit, estdat=siteQ)
  write.csv(inputMetrics, file.path(outputFolder, constitName, "inputs", paste0(constitSite, '.csv')), row.names=FALSE)
  
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
  allModels[[constitSite]] <- list(comp = comp, interpRect = interpRect, 
                                      rloadest5param = rloadest5param)
  
  #make predictions
  pred_rload <- predictSolute(rloadest5param, "flux", siteQ, 
                              se.pred = TRUE, date = TRUE)
  pred_interp <- predictSolute(interpRect, "flux", siteQ, 
                               se.pred = TRUE, date = TRUE)
  pred_comp <- predictSolute(comp, "flux", siteQ, se.pred = TRUE,
                             date = TRUE)
  nPreds <- nrow(siteQ)
  allPreds <- bind_rows(pred_rload, pred_interp, pred_comp)
  allPreds$model <- c(rep("rloadest",nPreds), rep("interp",nPreds), rep("composite", nPreds))
    
  #TODO: model metrics
  annualPreds <- bind_rows(
    summarizePreds(pred_rload, siteMeta, "total", model.name = "rloadest"),
    summarizePreds(pred_interp, siteMeta, "total", model.name = "interpolation"),
    summarizePreds(pred_comp, siteMeta, "total", model.name = "composite"))
  write.csv(x = annualPreds, file = file.path(outputFolder, constitName, "annual", paste0(constitSite, '.csv')), row.names=FALSE)
  
  #TODO: plots
  writePDFreport(file = file.path(outputFolder, constitName, paste(constitSite, "report.pdf", sep = "_")),
                 intdat = siteConstit, estdat = siteQ, allPreds = allPreds, 
                 meta = siteMeta, inputCSV = inputMetrics, annualCSV = annualPreds)
  
    
  #TODO: verbose option to print output?
  
  #TODO: compute and save whatever is supposed to go in the multiYear data.frame for this site-constituent combo
  # multiYearSummary <- ...
  # write.csv(x = multiYearSummary, file = file.path(outputFolder, constitName, "multiYear", paste0(constitSite, '.csv')), row.names=FALSE)
  
  message(paste('Finished processing constituent file', fileDF$constitFile[i], '\n'))
}

allInputs <- summarizeCsvs('inputs', fileDF, outputFolder) 
allAnnual <- summarizeCsvs('annual', fileDF, outputFolder) 
allMultiYear <- summarizeCsvs('multiYear', fileDF, outputFolder) 
