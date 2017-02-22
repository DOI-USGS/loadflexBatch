#' Convert EGRET files into the format currently being used by the Brazillian 
#' ANA for load analyses.
#' 
#' @param egret.files a character vector filenames of EGRET eList objects saved 
#'   in .RData format
#' @param ANA.dir the directory filepath where the converted files should be 
#'   written
#' @return TRUE if the conversion was successful. Throws an error otherwise.
EGRETtoANA <- function(egret.files, ANA.dir) {
  
  library(dplyr)
  library(EGRET)
  library(loadflex)
  
  # read the RData files and confirm that they're EGRET eLists
  elists <- lapply(egret.files, function(egret.file) {
    elist.name <- load(egret.file)
    egret.elist <- get(elist.name)
    if(class(egret.elist) != 'egret' ||
       !all.equal(names(egret.elist), c('INFO','Daily','Sample','surfaces'))) {
      stop(egret.file, " doesn't contain an eList with elements INFO, Daily, Sample, and surface")
    }
    return(egret.elist)
  })
  
  # extract the INFOs into a siteInfo table and write to csv
  constitInfo <- bind_rows(
    lapply(elists, function(elist) {
      elist$INFO %>%
        mutate(
          conc.units = {
            if(all(grepl('mg/l', tolower(param.units)))) {
              'mg L^-1'
            } else {
              stop('expected param.units of mg/l, per Table 2 of https://pubs.usgs.gov/tm/04/a10/pdf/tm4A10.pdf')
            }
          },
          matching.site = staAbbrev
        ) %>%
        select(
          matching.site,
          site.id = staAbbrev, # as in convertToEGRET
          site.name = shortName, # as in convertToEGRET
          lat = dec_lat_va,
          lon = dec_long_va,
          basin.area = drainSqKm,
          constituent = constitAbbrev, # as in convertToEGRET
          consti.name = paramShortName, # as in convertToEGRET
          units = conc.units # as in convertToEGRET. always mg/L
        ) %>% 
        return()
    })
  )
  # Discharge units: Table 1 of https://pubs.usgs.gov/tm/04/a10/pdf/tm4A10.pdf
  QInfo <- constitInfo %>%
    mutate(
      constituent = 'Q',
      consti.name = 'Discharge, cubic meters per second',
      units = 'cms'
    )
  siteInfo <- bind_rows(constitInfo, QInfo) %>% distinct()
  write.csv(siteInfo, file.path(ANA.dir, 'siteInfo.csv'), row.names=FALSE)
  
  # extract the constituent data and write to csvs
  lapply(elists, function(elist) {
    constituent <- elist$INFO$constitAbbrev
    const.dir <- file.path(ANA.dir, constituent)
    if(!dir.exists(const.dir)) dir.create(const.dir)
    site.id <- elist$INFO$staAbbrev
    constitData <- 
      elist$Sample %>%
      # next lines assume concentrations are always left censored or uncensored.
      # we know our first two files are entirely uncensored.
      mutate(
        constit = ifelse(Uncen == 1, ConcAve, ConcHigh),
        status = ifelse(Uncen == 1, 1, 2) # Uncen can be 1 (uncensored) or 0 (censored)
      ) %>% 
      select(date = Date, Q, constit, status) %>%
      setNames(c('date', 'Q', constituent, 'status'))
    write.csv(constitData, file.path(const.dir, paste0(site.id, '.csv')), row.names=FALSE)
  })
  
  # extract the discharge data and write to csvs
  lapply(elists, function(elist) {
    Q.dir <- file.path(ANA.dir, 'Q')
    if(!dir.exists(Q.dir)) dir.create(Q.dir)
    site.id <- elist$INFO$staAbbrev
    QData <- 
      elist$Daily %>%
      select(date = Date, Q)
    write.csv(QData, file.path(Q.dir, paste0(site.id, '.csv')), row.names=FALSE)
  })

  detach(package:loadflex, unload=TRUE) # because loadflex imports EGRET, and we need to detach EGRET
  detach(package:EGRET, unload=TRUE) # because getInfo masks loadflex::getInfo
  
  return(TRUE) # if we get this far, the conversion succeeded (otherwise we stop()ed above)
}

EGRETtoANA(egret.files=c('Hirsch_sites/EGRET/Musk.NO23.RData', 'Hirsch_sites/EGRET/Rac.NO23.RData'), ANA.dir='Hirsch_sites/input')
