#add censored status column to wq inputs
#need wd set where this script is
lf <- list.files(c('PT', "NO3"), full.names = TRUE)

library(dplyr)

for(file in lf) {
  df <- read.csv(file, stringsAsFactors = FALSE)  
  df <- mutate(df, status = rep(1, nrow(df)))
  
  #randomly make some 0 or 2
  df$status[sample(1:nrow(df), 5)] <- sample(c(0,2), 5, replace = TRUE)
  
  write.csv(file = file, x = df, row.names = FALSE)
}