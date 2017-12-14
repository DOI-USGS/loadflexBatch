library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)

#### NO3 for RONC02800 ####

predsNO3 <- bind_rows(
  # with the original DECTIME and DECTIME2
  read_csv('three_ANA_sites/output/NO3/monthly/RONC02800.csv') %>%
    mutate(fixed = 'none'),
  # with fixed DECTIME and DECTIME2 to various base years
  read_csv('three_ANA_sites/output_f2012/NO3/monthly/RONC02800.csv') %>%
    mutate(fixed = '2012'),
  read_csv('three_ANA_sites/output_f2007/NO3/monthly/RONC02800.csv') %>%
    mutate(fixed = '2007')
)
inputNO3 <- read_csv('three_ANA_sites/input/NO3/RONC02800.csv') %>%
  mutate(flux = NO3 * Q * loadflex::flowconcToFluxConversion('m^3 s^-1', 'mg L^-1', 'kg y^-1'))

ggplot(predsNO3, aes(x=as.Date(sprintf("%s-01", Month)), y=RL7.Flux_Rate)) +
  geom_step(aes(color=fixed, alpha='linealpha')) +
  geom_point(data=inputNO3, aes(x=as.Date(date), y=flux), color='gold3', alpha=0.7) +
  ggtitle("NO3, RONC02800, RL7") + xlab("Date") + ylab("Flux Rate (kg y^-1)") +
  scale_color_discrete("Base Year") +
  scale_alpha_manual(values=c(linealpha=0.7), guide=FALSE) +
  theme_bw()
ggsave('example_code/fixed_month_NO3_RONC02800_RL7.png', width=9, height=4)

ggplot(predsNO3, aes(x=as.Date(sprintf("%s-01", Month)), y=RL5.Flux_Rate)) +
  geom_step(aes(color=fixed, alpha='linealpha')) +
  geom_point(data=inputNO3, aes(x=as.Date(date), y=flux), color='gold3', alpha=0.7) +
  ggtitle("NO3, RONC02800, RL5") + xlab("Date") + ylab("Flux Rate (kg y^-1)") +
  scale_color_discrete("Base Year") +
  scale_alpha_manual(values=c(linealpha=0.7), guide=FALSE) +
  theme_bw()
ggsave('example_code/fixed_month_NO3_RONC02800_RL5.png', width=9, height=4)


#### PT for RONC02800 ####

predsPT <- bind_rows(
  # with the original DECTIME and DECTIME2
  read_csv('three_ANA_sites/output/PT/monthly/RONC02800.csv') %>%
    mutate(fixed = 'none'),
  # with fixed DECTIME and DECTIME2 to various base years
  read_csv('three_ANA_sites/output_f2012/PT/monthly/RONC02800.csv') %>%
    mutate(fixed = '2012'),
  read_csv('three_ANA_sites/output_f2007/PT/monthly/RONC02800.csv') %>%
    mutate(fixed = '2007')
)
inputPT <- read_csv('three_ANA_sites/input/PT/RONC02800.csv') %>%
  mutate(flux = PT * Q * loadflex::flowconcToFluxConversion('m^3 s^-1', 'mg L^-1', 'kg y^-1'))

ggplot(predsPT, aes(x=as.Date(sprintf("%s-01", Month)), y=RL5.Flux_Rate)) +
  geom_step(aes(color=fixed, alpha='linealpha')) +
  geom_point(data=inputPT, aes(x=as.Date(date), y=flux), color='gold3', alpha=0.7) +
  ggtitle("PT, RONC02800, RL5") + xlab("Date") + ylab("Flux Rate (kg y^-1)") +
  scale_color_discrete("Base Year") +
  scale_alpha_manual(values=c(linealpha=0.7), guide=FALSE) +
  theme_bw()
ggsave('example_code/fixed_month_PT_RONC02800_RL5.png', width=9, height=4)

ggplot(predsPT, aes(x=as.Date(sprintf("%s-01", Month)), y=RL7.Flux_Rate)) +
  geom_step(aes(color=fixed, alpha='linealpha')) +
  geom_point(data=inputPT, aes(x=as.Date(date), y=flux), color='gold3', alpha=0.7) +
  ggtitle("PT, RONC02800, RL7") + xlab("Date") + ylab("Flux Rate (kg y^-1)") +
  scale_color_discrete("Base Year") +
  scale_alpha_manual(values=c(linealpha=0.7), guide=FALSE) +
  theme_bw()
ggsave('example_code/fixed_month_PT_RONC02800_RL7.png', width=9, height=4)

#### multiYear preds - NO3 ####

MYPredsNO3 <- bind_rows(
  # with the original DECTIME and DECTIME2
  read_csv('three_ANA_sites/output/NO3/multiYear/RONC02800.csv') %>%
    mutate(fixed = 'none'),
  # with fixed DECTIME and DECTIME2 to various base years
  read_csv('three_ANA_sites/output_f2012/NO3/multiYear/RONC02800.csv') %>%
    mutate(fixed = '2012'),
  read_csv('three_ANA_sites/output_f2007/NO3/multiYear/RONC02800.csv') %>%
    mutate(fixed = '2007')
)

dat <- full_join(
  gather(MYPredsNO3, MFR, Flux_Rate, RL5.Flux_Rate, RL7.Flux_Rate) %>%
    mutate(Model=sapply(strsplit(MFR, "\\."), `[[`, 1)) %>%
    select(site.id, Model, fixed, Flux_Rate),
  full_join(
    gather(MYPredsNO3, MFR, CI_lower, RL5.CI_lower, RL7.CI_lower) %>%
      mutate(Model=sapply(strsplit(MFR, "\\."), `[[`, 1)) %>%
      select(site.id, Model, fixed, CI_lower),
    gather(MYPredsNO3, MFR, CI_upper, RL5.CI_upper, RL7.CI_upper) %>%
      mutate(Model=sapply(strsplit(MFR, "\\."), `[[`, 1)) %>%
      select(site.id, Model, fixed, CI_upper),
    by=c('site.id','Model','fixed')),
  by=c('site.id','Model','fixed'))
dodge <- position_dodge(width=0.9)
ggplot(dat, aes(x=Model, y=Flux_Rate, fill=fixed, group=fixed)) +
  geom_col(aes(alpha='baralpha'), position=dodge) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper, color=fixed), position=dodge, width=0.4) +
  ggtitle("NO3, RONC02800, 2 regression models and 3 base years") + 
  xlab("Model") + ylab("Flux Rate (kg y^-1)") +
  scale_color_discrete(guide=FALSE) +
  scale_alpha_manual(values=c(baralpha=0.4), guide=FALSE) +
  theme_bw()
ggsave('example_code/fixed_multiYear_NO3_RONC02800.png', width=8, height=4)

#### multiYear preds - PT ####

MYPredsPT <- bind_rows(
  # with the original DECTIME and DECTIME2
  read_csv('three_ANA_sites/output/PT/multiYear/RONC02800.csv') %>%
    mutate(fixed = 'none'),
  # with fixed DECTIME and DECTIME2 to various base years
  read_csv('three_ANA_sites/output_f2012/PT/multiYear/RONC02800.csv') %>%
    mutate(fixed = '2012'),
  read_csv('three_ANA_sites/output_f2007/PT/multiYear/RONC02800.csv') %>%
    mutate(fixed = '2007')
)

dat <- full_join(
  full_join(
    gather(MYPredsPT, MFR, Flux_Rate, RL5.Flux_Rate, RL7.Flux_Rate) %>%
      mutate(Model=sapply(strsplit(MFR, "\\."), `[[`, 1)) %>%
      select(site.id, Model, fixed, Flux_Rate),
    gather(MYPredsPT, MFR, SE, RL5.SE, RL7.SE) %>%
      mutate(Model=sapply(strsplit(MFR, "\\."), `[[`, 1)) %>%
      select(site.id, Model, fixed, SE),
    by=c('site.id','Model','fixed')),
  full_join(
    gather(MYPredsPT, MFR, CI_lower, RL5.CI_lower, RL7.CI_lower) %>%
      mutate(Model=sapply(strsplit(MFR, "\\."), `[[`, 1)) %>%
      select(site.id, Model, fixed, CI_lower),
    gather(MYPredsPT, MFR, CI_upper, RL5.CI_upper, RL7.CI_upper) %>%
      mutate(Model=sapply(strsplit(MFR, "\\."), `[[`, 1)) %>%
      select(site.id, Model, fixed, CI_upper),
    by=c('site.id','Model','fixed')),
  by=c('site.id','Model','fixed')) %>%
  mutate(
    CI_lower2 = Flux_Rate - 1.96*SE,
    CI_upper2 = Flux_Rate + 1.96*SE)
dodge <- position_dodge(width=0.9)
ggplot(dat, aes(x=Model, y=Flux_Rate, fill=fixed, group=fixed)) +
  geom_col(aes(alpha='baralpha'), position=dodge) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper, color=fixed), position=dodge, width=0.4) +
  # geom_errorbar(aes(ymin=CI_lower2, ymax=CI_upper2, color=fixed), position=dodge, width=0.3, linetype='dashed') +
  ggtitle("PT, RONC02800, 2 regression models and 3 base years") + 
  xlab("Model") + ylab("Flux Rate (kg y^-1)") +
  scale_color_discrete(guide=FALSE) +
  scale_alpha_manual(values=c(baralpha=0.4), guide=FALSE) +
  theme_bw()
ggsave('example_code/fixed_multiYear_PT_RONC02800.png', width=8, height=4)

