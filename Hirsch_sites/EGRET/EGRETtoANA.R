library(dplyr)
library(EGRET)

EGRETtoANA <- function(egret.files, ANA.dir) {
  
  # read the RData files and confirm that they're EGRET eLists
  elists <- lapply(egret.files, function(egret.file) {
    elist.name <- load(egret.file)
    egret.elist <- get(elist.name)
    if(!all.equal(names(egret.elist), c('INFO','Daily','Sample','surfaces')) || class(egret.elist) != 'egret') {
      stop("are you sure this file contains an eList?")
    }
    return(egret.elist)
  })
  
  # extract the INFOs into a siteInfo table and write to csv
  constitInfo <- bind_rows(
    lapply(elists, function(elist) {
      elist$INFO %>%
        mutate(
          conc.units = c("mg/l as N"='mg L^-1')[param.units],
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
          units = conc.units # as in convertToEGRET
        )
    })
  )
  QInfo <- constitInfo %>%
    mutate(
      constituent = 'Q',
      consti.name = 'Discharge, cubic feet per second', # assuming parameter_cd 00060
      units = 'cfs'
    )
  siteInfo <- bind_rows(constitInfo, QInfo)
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
        status = ifelse(Uncen == 1, 1, 2)
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
}

EGRETtoANA(c('Hirsch_sites/EGRET/Musk.NO23.RData', 'Hirsch_sites/EGRET/Rac.NO23.RData'), ANA.dir='Hirsch_sites/input')

detach(package:loadflex, unload=TRUE) # because loadflex imports EGRET, and we need to detach EGRET
detach(package:EGRET, unload=TRUE) # because getInfo masks loadflex::getInfo

