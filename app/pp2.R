# Preprocessing of the raw data
# launched with 
# R CMD BATCH --vanilla --no-restore app/preprocessing.R app/log_pp.txt

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
                        week = lubridate::week(date))

# stations we look at
stns_id <- c(568, 4333, 5051)
p <- 3
data <- data %>%
  filter(station_id %in% stns_id) %>%
  arrange(station_id, date)

# quick look at summary of the data (clearly shows time trends)
data0 <- data %>% group_by(station_name, year) %>%
  summarise(mm_temp = mean(mean_temp, na.rm=T))
ggplot(data0, aes(x=year, y=mm_temp)) +
  geom_line() +
  geom_point() +
  facet_wrap(~station_name)


# Estimate mean conditional on time of year -------------------------------
year_len <- 365 + 6/24 + 9/60/24 + 9/60^2/24

# create sines-cosines basis (for seasonal trend)
nSC <- 12
s_cols <- paste0("s", 1:nSC)
c_cols <- paste0("c", 1:nSC)
for (i in 1:nSC) data <- data %>%
  mutate(!!s_cols[i] := sin(i*2*pi*time/year_len),
         !!c_cols[i] := cos(i*2*pi*time/year_len))


# create natural spline basis (for long term time trend)
# for this check range for each station
row_ids <- sapply(stns_id, \(id) data$station_id == id)
bss <- lapply(1:p, \(k){
  data[row_ids[,k],] %>% filter(day == 21, month %in% c(6,12)) %>% select(time) %>% unlist
})

for(b in paste("b", 1:max(sapply(bss, length)))) data <- data %>% mutate(!!b := NA)
nB <- rep(0, p); b_cols <- list()
for(k in 1:p){
  B <- splines::bs(data[row_ids[,k],]$time, knots = c(bss[[k]]), degree = 2)
  B <- B[,colSums(abs(B[!is.na(data[row_ids[,k],]$mean_temp),])) > 1e-5]
  nB[k] <- ncol(B)
  b_cols[[length(b_cols)+1]] <- paste0("b", 1:nB[k])
  data[row_ids[,k],b_cols[[k]]] <- B
  rm(B)
}


# fit least-squares
non_na <- !is.na(data$mean_temp) # for later
lms <- lapply(1:p, \(k){
  x_cols <- c("mean_temp", s_cols, c_cols, b_cols[[k]])
  lm(formula = "mean_temp~.", data = data[row_ids[,k],x_cols])
})

# Bakerville
summary(lms[[1]])
# Ottawa
summary(lms[[2]])
# Toronto
summary(lms[[3]])

# register the results
data <- data %>% mutate(mu = NA , ctemp = NA)
for(k in 1:p){
  data[non_na & row_ids[,k],] <- data[non_na & row_ids[,k],] %>% 
    mutate(mu = lms[[k]]$fitted.values, ctemp = mean_temp - mu)
}

# remove b cols of data (for memory)
data <- data %>% 
  select(!any_of(b_cols[[which.max(sapply(b_cols, length))]])) %>%
  select(!any_of(c(s_cols, c_cols)))
gc()

# quick check
data0 <- data %>% group_by(year, station_name) %>%
  summarise(m_mt = mean(mean_temp, na.rm=T),
            m_mu = mean(mu, na.rm=T))
ggplot(data0, aes(x=year, y=m_mt, col=station_name, group=station_name)) +
  theme_light() +
  geom_line(data = data0, aes(y=m_mu), col=1) +
  geom_point()

data0 <- data %>% filter(year %in% 1950:1970) %>%
  group_by(year, month, station_name) %>%
  summarise(m_mt = mean(mean_temp, na.rm=T),
            m_mu = mean(mu, na.rm=T))
ggplot(data0, aes(x=year+month/12, y=m_mt, col=station_name, group=station_name)) +
  geom_line() + geom_point() +
  geom_line(data = data0, aes(y=m_mu)) +
  facet_wrap(~station_name)


# random check of the fits
# some plots might show a break (when the data goes from Ottawa to Toronto)
k <- sample(nrow(data),1)
ii <- k + 1:365
plot(data$time[ii], data$mean_temp[ii])
lines(data$time[ii], data$mu[ii], col=2)


# note that this breaks the equalities in the data
# but there is still some heterogeneity in the variance
N <- 20000 
par(mfrow=c(p,2), mar=c(1,1,0,0))
for(k in 1:p){
  plot(data$yday[row_ids[,k]][1:N], data$mean_temp[row_ids[,k]][1:N], cex=.25)
  plot(data$yday[row_ids[,k]][1:N], data$ctemp[row_ids[,k]][1:N], cex=.25)
}



# Compute pseudo-observations using ECDF ----------------------------------
data <- data %>% mutate(pseudo_temp = NA, diff = NA, w1 = NA, w = NA)
for(k in 1:p){
  id <- stns_id[k]
  cat("Working on ECDF for station ", id, ".\n")
  row_id <- row_ids[,k]
  
  for(tt in 1:366){
    if(tt %% 20 == 0) cat("progess: ", round(tt/366*100,2), "%\n")
    
    # construct weights based on yday using gaussian kernel
    data[row_id,] <- data[row_id,] %>%
      mutate(diff = tt - (time %% year_len)) %>%
      mutate(diff = abs(pmin(diff, year_len - diff))) %>%
      mutate(w1 = dnorm(diff,0,2)) %>%
      mutate(w = w1/sum(w1, na.rm = T)) 
    
    # construct ecdf function using all the data
    wecdf <- spatstat.geom::ewcdf(data[row_id,]$ctemp, data[row_id,]$w)
    
    # compute ecdf only for yday we care for
    ind <- which(data[row_id,]$yday == tt & !is.na(data[row_id,]$ctemp))
    data[row_id,][ind,] <- data[row_id,][ind,] %>%
      mutate(pseudo_temp = wecdf(ctemp))
    
    # heuristic correction for the largest (pseudo_obs = 1)
    # to avoid modifying all values (say, with pseudo_obs = pseudo_obs*n_ind/(n_ind+1))
    n_ind <- length(ind)
    data[row_id,][ind,] <- data[row_id,][ind,] %>%
      mutate(pseudo_temp = ifelse(pseudo_temp == 1, 1-1/(2*n_ind), pseudo_temp))
  }
}



# resulting histograms are even better than expected
ggplot(data, aes(x = pseudo_temp, fill=as.factor(yday))) +
  theme(legend.position = "none") +
  geom_histogram(breaks=seq(0,1,.05), position = position_stack()) +
  facet_wrap(~station_name)

# quick check of the other operations (tt = 366, last it. of for loop)
s <- sample(nrow(data) - 1000,1)
par(mfrow=c(3,1))
plot(data$pseudo_temp[s + 1:100], type="o") # makes sense
plot(data$date[s + 1:1000], data$diff[s + 1:1000]) # linear in time, as expected
plot(data$date[s + 1:1000], data$w[s + 1:1000], type="l") # select only end of year periods

# note that there are equalities, but very very few of them. can disregard.
apply(row_ids, 2, \(r){
  tab <- table(data[r,]$pseudo_temp)
  tab[tab > 1]
})



# Save for later use ------------------------------------------------------
saveRDS(data, "app/data/data_pp2.rds")
