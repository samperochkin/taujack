# Data analysis


# packages ----------------------------------------------------------------
library(tidyverse)



# load data ---------------------------------------------------------------
stns_id <- readRDS("app/data/stns_id.rds")
data <- readRDS("app/data/data_res.rds")


# quick analysis of univariate residuals ----------------------------------

# uniformly distributed?
hist(data$res, probability = T, breaks=50)

par(mfrow = c(2,4), mar = c(2,2,3,1))
for(stn in stns_id){
  ind <- which(data$station_id == stn)
  hist(data[ind,]$res, main = paste0("stn: ", stn),
       probability = T, breaks=50)
}

# make sure plot area is large enough
# seems okay, some extreme values, but overall not too bad
par(mfrow = c(2,3), mar = c(2,2,3,1))
for(stn in stns_id){
  ind <- which(data$station_id == stn)
  hist(data[ind,]$res, main = paste0("stn: ", stn),
         probability = T, breaks=30)
}

# acf and pacf quick check
# done on about 9 years of data (uninterupted, i.e., no NAs)
# conclusion: this seems satisfactory for illustrative purposes

# check for various lags
lag_max = 60; ylim = c(-.1,.1)
# lag_max = 24*30; ylim = c(-.05,.05)
# lag_max = 24*365*2; ylim = c(-.025,.025)

par(mfrow = c(2,3), mar = c(2,2,3,1))
for(stn in stns_id){
  ind <- which(data$station_id == stn & !is.na(data$res)) # this is where we subset
  d <- diff(ind)
  ks <- which(d > 1)
  dd <- diff(ks)
  m <- ks[which.max(dd)-1]
  k <- ks[m]
  
  ind <-  ind[k + 1:dd[m]]
  table(is.na(data[ind,]$res))
  
  
  k <- 1
  N <- 2000
  while(any(diff(ind[(k-1)*N + 1:5000]) > 1)) k <- k+1
  while(any(diff(ind[(k-1)*N + 1:5000]) > 1)) k <- k+1
  acf(qnorm(data[ind[(k-1)*N + 1:5000],]$res), lag.max = lag_max, ylim = ylim,
      main = paste0("stn: ", stn))
}

for(stn in stns_id){
  ind <- which(data$station_id == stn & !is.na(data$temp)) # this is where we subset
  k <- 1
  while(any(diff(ind[(k-1)*2000 + 1:2000]) > 1)) k <- k+1
  pacf(qnorm(data[ind[(k-1)*2000 + 1:2000],]$res), lag.max = lag_max, ylim = ylim,
       main = paste0("stn: ", stn))
}

