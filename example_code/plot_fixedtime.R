monthlySummary # with the original DECTIME and DECTIME2
monthlySummaryF # with fixed DECTIME and DECTIME2
siteConstit # NO3 for RONC02800 shows pretty impressive effects

ggplot(monthlySummary, aes(x=as.Date(sprintf("%s-01", Month)), y=RL7.Flux_Rate)) +
  geom_line() +
  geom_line(color='blue', data=monthlySummaryF) +
  geom_point(data=siteConstit, aes(x=as.Date(date), y=NO3*Q*365))

