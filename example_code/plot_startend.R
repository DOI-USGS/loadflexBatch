library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)

#### NO3 for RONC02800, clipped at 10/1/2002 or 10/1/2012 ####

predsNO3 <- bind_rows(
  read_csv('three_ANA_sites/output/NO3/annual/RONC02800.csv') %>%
    mutate(date.start = 'NA'),
  read_csv('three_ANA_sites/output_RONCstart/NO3/annual/RONC02800.csv') %>%
    mutate(date.start = '2012-10-01')
) %>%
  mutate(date.start = ordered(date.start, levels=c('NA','2012-10-01')))

ggplot(predsNO3, aes(x=Water_Year, y=RL5.Flux_Rate, fill=date.start, linetype=date.start)) +
  geom_col(position='dodge') +
  xlab("Water Year") + ylab("Flux Rate (kg y^-1) from RL5") + ggtitle("NO3, RONC028000, annual fluxes") +
  scale_fill_manual(values=c('NA'='orange','2012-10-01'='seagreen')) +
  theme_bw()
ggsave('example_code/startend_NO3_RONC02800_RL5.png', width=7, height=4)
