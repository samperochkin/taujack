
# packages ----------------------------------------------------------------
library(tidyverse)
library(gamlss)
library(parallel)



# load and preprocess data ------------------------------------------------
filepaths <- list.files("data", pattern = "temp_data", full.names = T)
stns_id <- readRDS("data/stns_id.rds")

filepaths <- filepaths[grepl(stns_id[1], filepaths)]


data <- lapply(filepaths, readRDS) |> dplyr::bind_rows()
data <- data %>% mutate(date = as.Date(date),
                        year = as.integer(year),
                        month = as.integer(month),
                        day = as.integer(day),
                        hour = lubridate::hour(strptime(hour, format = '%H:%M')))
data <- data %>% mutate(time = as.integer(data$date) + data$hour/24,
                        yday = lubridate::yday(date),
                        week = lubridate::week(date),
                        seas = lunar::terrestrial.season(date))
str(data)



# setup for gamlss model --------------------------------------------------
table(round(diff(data$time),10))
n <- nrow(data)

# create lagged variables
data$temp_lag1 <- c(rep(NA,1), data$temp[-n])
data$temp_lag2 <- c(rep(NA,2), data$temp[-((n-1):n)])
data$temp_lag3 <- c(rep(NA,3), data$temp[-((n-2):n)])
data$temp_lag6 <- c(rep(NA,6), data$temp[-((n-5):n)])
data$temp_lag12 <- c(rep(NA,12), data$temp[-((n-11):n)])
data$temp_lag24 <- c(rep(NA,24), data$temp[-((n-23):n)])

# to identify non-NA rows later
data <- data %>% mutate(keep = !is.na(temp) & 
                          !is.na(temp_lag1) &
                          !is.na(temp_lag2) &
                          !is.na(temp_lag3) &
                          !is.na(temp_lag6) &
                          !is.na(temp_lag12) &
                          !is.na(temp_lag24))

# create time variables
data <- data[order(data$station_id),]
year_len <- 365 + 6/24 + 9/60/24 + 9/60^2/24
moon_len <- 1 + 0/24 + 50/60/24

data <- data %>% mutate(s1_time = sin(2*pi*time/year_len),
                        c1_time = cos(2*pi*time/year_len),
                        s2_time = sin(4*pi*time/year_len),
                        c2_time = cos(4*pi*time/year_len),
                        s3_time = sin(6*pi*time/year_len),
                        c3_time = cos(6*pi*time/year_len))

# subset data so that no NA creeps up in usused columns
data0 <- data[, c("keep", "temp", "date", "time", "month", "day", "hour",
                  grep("temp_lag", names(data), value = T),
                  grep("_time", names(data), value = T))]
ind <- which(data0$keep)

ff <- temp ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12 + temp_lag24)*
  (s1_time + c1_time + s2_time + c2_time + s3_time + c3_time)

mclapply(0:23, \(h){
  indh <- ind[data0[ind,]$hour == h]
  mod <- gamlss(formula = ff, ff, ff, ff,
                family = gamlss.dist::ST2(),
                data = data0[indh, 1:19],
                control = gamlss.control(c.crit = .01, n.cyc = 200),
                i.control = glim.control(cc = .01))
  
  saveRDS(cbind(mod$mu.coefficients,
                mod$sigma.coefficients,
                mod$nu.coefficients,
                mod$tau.coefficients), paste0("app/P_", h, ".rds"))
}, mc.cores = 12)
