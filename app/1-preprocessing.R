# Preprocessing of the raw data
# launched with 
# R CMD BATCH --no-restore app/preprocessing.R app/log_pp.txt

# packages ----------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(lunar)
library(ggplot2)
library(spatstat.geom) # weighted ecdf



# load and preprocess data ------------------------------------------------
stns_id <- readRDS("app/data/stns_id.rds")
weathercan::stations() %>% filter(station_id %in% stns_id, interval == "day") %>%
  print(nrow = nrow(.))

filepaths <- paste0("app/data/daily_data_", stns_id, ".rds")
data <- lapply(filepaths, readRDS) |> dplyr::bind_rows()
str(data)

# time variables setup
data <- data %>% mutate(date = as.Date(date),
                        year = as.integer(year),
                        month = as.integer(month),
                        day = as.integer(day))

data <- data %>% mutate(time = as.integer(date),
                        yday = lubridate::yday(date),
                        week = lubridate::week(date),
                        seas = lunar::terrestrial.season(date))

# stations we look at
stns_id <- c(5051)
data <- data %>%
  filter(station_id %in% stns_id) %>%
  arrange(station_id, date)




# Estimate mean conditional on time of year -------------------------------
year_len <- 365 + 6/24 + 9/60/24 + 9/60^2/24

# create sines-cosines basis
nSC <- 8
for (i in 1:nSC){
  s <- paste0("s", i)
  c <- paste0("c", i)
  data <- data %>%
    mutate(!!s := sin(i*2*pi*time/year_len),
           !!c := cos(i*2*pi*time/year_len))
}

# create natural spline basis for time trend
nB <- 8
B <- splines::ns(data$time, df = nB)
for (i in 1:nB){
  b <- paste0("b", i)
  data <- data %>%
    mutate(!!b := B[,i])
}

# fit least-squares
non_na <- which(!is.na(data$max_temp))
ff <- paste0(
  "max_temp ~ b1", paste0(" + b", 1:nB, collapse = ""),
  paste0(" + s", 1:nSC, collapse = ""),
  paste0(" + c", 1:nSC, collapse = "")
) %>% as.formula

lm_fit <- lm(formula = ff, data = data[non_na,])
summary(lm_fit)

data <- data %>% mutate(mu = NA , ctemp = NA)
data[non_na,] <- data[non_na,] %>% 
  mutate(mu = lm_fit$fitted.values, ctemp = max_temp - mu)

data0 <- data %>% group_by(year) %>%
  summarise(m_mt = mean(max_temp, na.rm=T),
            m_mu = mean(mu, na.rm=T))
plot(data0$year, data0$m_mt, type="o")  
lines(data0$year, data0$m_mu, col=2)  

data0 <- data %>% filter(year %in% 1950:1970) %>%
  group_by(year, month) %>%
  summarise(m_mt = mean(max_temp, na.rm=T),
            m_mu = mean(mu, na.rm=T))
plot(data0$year+data0$month/12, data0$m_mt, type="o")  
lines(data0$year+data0$month/12, data0$m_mu, col=2)  


k <- sample(nrow(data),1)
ii <- k + 1:365
plot(data$time[ii], data$max_temp[ii])
lines(data$time[ii], data$mu[ii], col=2)


# note that this breaks the equalities in the data
# but there is still some heterogeneity in the variance
plot(data$yday[1:20000], data$max_temp[1:20000], cex=.25)
plot(data$yday[1:20000], data$ctemp[1:20000], cex=.25)



# Compute pseudo-observations using ECDF ----------------------------------
data <- data %>% mutate(pseudo_temp = NA)
for(tt in 1:366){
  
  if(tt %% 20 == 0) cat("progess: ", round(tt/366*100,2), "%\n")
  
  # construct weights based on yday using gaussian kernel
  data <- data %>%
    mutate(diff = tt - (time %% year_len)) %>%
    mutate(diff = abs(pmin(diff, year_len - diff))) %>%
    mutate(w1 = dnorm(diff,0,1)) %>%
    mutate(w = w1/sum(w1, na.rm = T)) 
  
  # construct ecdf function using all the data
  wecdf <- spatstat.geom::ewcdf(data$ctemp, data$w)
  
  # compute ecdf only for yday we care for the appropriate subset
  ind <- which(data$yday == tt & !is.na(data$ctemp))
  data[ind,] <- data[ind,] %>%
    mutate(pseudo_temp = wecdf(ctemp))
}

# resulting histogram is even better than expected
ggplot(data, aes(x = pseudo_temp, fill=as.factor(yday))) +
  theme(legend.position = "none") +
  geom_histogram(breaks=seq(0,1,.05), position = position_stack())

# quick check of the other operations (tt = 366, last it. of for loop)
s <- sample(nrow(data) - 1000,1)
plot(data$pseudo_temp[s + 1:100], type="o") # makes sense
plot(data$date[s + 1:1000], data$diff[s + 1:1000]) # linear in time, as expected
plot(data$date[s + 1:1000], data$w[s + 1:1000], type="l") # select only end of year periods

# note that there are equalities, but very very few of them. can disregard.
tab <- table(data$pseudo_temp)
tab[tab > 1]
sum(tab)



# Save for later use ------------------------------------------------------
saveRDS(data, "app/data/data_pp.rds")
