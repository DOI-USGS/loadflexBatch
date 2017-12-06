# Example code for manually producing time series plots of predictions at
# various resolutions

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(lubridate)

control_file <- 'three_ANA_sites.yml'
inputs <- readInputs(control_file)
constitName <- 'NO3'
sites <- c("RONC02800") # "MOGU02900", "ORIZ02900"

outdir <- file.path(inputs$outputFolder, constitName)
preds <- lapply(file.path(outdir, intersect(c('monthly','seasonal','annual','multiYear'), dir(outdir))) %>% setNames(., nm=basename(.)), function(preddir) {
  lapply(dir(preddir, full.names=TRUE) %>% setNames(., nm=tools::file_path_sans_ext(basename(.))), function(predfile) {
    message(predfile)
    readr::read_csv(predfile) %>%
      {gather(., 'variable', 'value', 4:(ncol(.)-2))} %>%
      mutate(
        model = sapply(strsplit(variable, '\\.'), `[`, 1),
        var = sapply(strsplit(variable, '\\.'), function(var) paste(var[2:length(var)], collapse='.')),
        value = as.numeric(value)) %>%
      select(-variable) %>%
      spread(var, value)
  })
})

# prepare predictions
monthly <- bind_rows(preds$monthly) %>%
  mutate(Month = as.Date(sprintf('%s-01', Month))) %>%
  filter(site.id %in% sites) %>%
  arrange(constituent, site.id, Month)
seasonal <- bind_rows(preds$seasonal) %>%
  mutate(Season = as.Date(Season)) %>%
  filter(site.id %in% sites) %>%
  arrange(constituent, site.id, Season)
annual <- bind_rows(preds$annual) %>%
  mutate(Water_Year = as.Date(sprintf('%s-10-01', Water_Year-1))) %>%
  filter(site.id %in% sites) %>%
  arrange(constituent, site.id, Water_Year)

dir.create('temp', showWarnings = FALSE)

# plot monthly predictions
ggplot(monthly, aes(x=Month, y=Flux_Rate, color=model)) +
  geom_step() + #geom_point() +
  facet_wrap(~site.id, scales='free') + theme_bw()
ggsave('temp/monthly.png', width=10, height=4)

# plot seasonal predictions
ggplot(seasonal, aes(x=Season, y=Flux_Rate, color=model)) +
  geom_step() + #geom_point(aes(x=Season + as.difftime(45, units='days'))) +
  facet_wrap(~site.id, scales='free') + theme_bw()
ggsave('temp/seasonal.png', width=10, height=4)

# plot annual predictions
ggplot(annual, aes(x=Water_Year, y=Flux_Rate, color=model)) +
  geom_step() +# geom_point(aes(x=Water_Year + as.difftime(183, units='days'))) +
  facet_wrap(~site.id, scales='free') + theme_bw()
ggsave('temp/annual.png', width=10, height=4)

# plot annual + monthly predictions
ggplot(monthly, aes(x=Month, y=Flux_Rate, color=model)) +
  geom_step() + 
  geom_step(data=annual, aes(x=Water_Year)) +
  facet_wrap(~site.id, scales='free') + theme_bw()
ggsave('temp/monthlyannual.png', width=10, height=4)

# plot seasonal + monthly predictions
ggplot(filter(monthly, Month > as.Date('2011-09-30')), aes(x=Month, y=Flux_Rate, color=model)) +
  geom_step(data=filter(seasonal, Season > as.Date('2011-09-30')), aes(x=Season)) +
  geom_line() + 
  facet_wrap(~site.id, scales='free') + theme_bw()
ggsave('temp/monthlyseasonal.png', width=7, height=4)
