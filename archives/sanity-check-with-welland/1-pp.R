# Preprocessing of the raw data
# launched with 
# R CMD BATCH --vanilla --no-restore app/preprocessing.R app/log_pp.txt

# packages ----------------------------------------------------------------
library(tidyverse)
library(lubridate)
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
# stns_id <- c(4712, 5051)
stns_id <- c(4859, 4862)
data <- data %>%
  filter(station_id %in% stns_id) %>%
  arrange(station_id, date)

# quick look at summary of the data (clearly shows time trends)
data0 <- data %>% group_by(station_name, year) %>%
  summarise(mm_temp = mean(max_temp, na.rm=T))
ggplot(data0, aes(x=year, y=mm_temp)) +
  geom_line() +
  geom_point() +
  facet_wrap(~station_name)


# Estimate mean conditional on time of year -------------------------------
year_len <- 365 + 6/24 + 9/60/24 + 9/60^2/24

# create sines-cosines basis (for seasonal trend)
nSC <- 8
for (i in 1:nSC){
  s <- paste0("s", i)
  c <- paste0("c", i)
  data <- data %>%
    mutate(!!s := sin(i*2*pi*time/year_len),
           !!c := cos(i*2*pi*time/year_len))
}

# create natural spline basis (for long term time trend)
# for this check range for each station
row_id1 <- data$station_id == stns_id[1]
row_id2 <- data$station_id == stns_id[2]
r1 <- data[!is.na(data$max_temp) & row_id1,]$time %>% range %>% diff
r2 <- data[!is.na(data$max_temp) & row_id2,]$time %>% range %>% diff
# let us use 8 and 10 df, respectively, then.
nB <- round(c(r1,r2)/year_len*12) # ************************** RE-SPECIFY FOR EACH NEW STNS
B1 <- splines::ns(data[row_id1,]$time, df = nB[1])
B2 <- splines::ns(data[row_id2,]$time, df = nB[2])
colnames(B1) <- paste0("b", 1:nB[1])
colnames(B2) <- paste0("b", 1:nB[2])

data1 <- merge.data.frame(data[row_id1,], data.frame(B1, date=data[row_id1,]$date), by="date")
data2 <- merge.data.frame(data[row_id2,], data.frame(B2, date=data[row_id2,]$date), by="date")
# data <- rbind(data1,data2)
# 
# for (i in 1:max(ncol(B1), ncol(B2))){
#   if(i %% 10 == 0) cat("Progress:", round(100*i/max(ncol(B1), ncol(B2)),2), "%\n")
#   b <- paste0("b", i)
#   data <- data %>% mutate(!!b := NA)
#   if(i <= ncol(B1)) data[row_id1,] <- data[row_id1,] %>% mutate(!!b := B1[,i])
#   if(i <= ncol(B2)) data[row_id2,] <- data[row_id2,] %>% mutate(!!b := B2[,i])
# }

# fit least-squares
non_na <- !is.na(data$max_temp)
ss <- paste0(" + s", 1:nSC, collapse = "")
cs <- paste0(" + c", 1:nSC, collapse = "")
ff1 <- paste0("max_temp ~ b1", paste0(" + b", 2:ncol(B1), collapse = ""), ss, cs, collapse = "")
ff2 <- paste0("max_temp ~ b1", paste0(" + b", 2:ncol(B2), collapse = ""), ss, cs, collapse = "")
ff1 <- as.formula(ff1)
ff2 <- as.formula(ff2)

cols1 <- c("max_temp", paste0("s", 1:nSC), paste0("c", 1:nSC), paste0("b", 1:nB[1]))
cols2 <- c("max_temp", paste0("s", 1:nSC), paste0("c", 1:nSC), paste0("b", 1:nB[2]))
lm_fit1 <- lm(formula = ff1, data = data1[,cols1])
lm_fit2 <- lm(formula = ff2, data = data2[,cols2])
summary(lm_fit1)
summary(lm_fit2)
# # stns 1
# lm_fit1 <- lm(formula = ff1, data = data[non_na,] %>% filter(station_id == stns_id[1]))
# lm_fit1 <- lm(formula = ff1, data = data[non_na,] %>% filter(station_id == stns_id[1]))
# summary(lm_fit1)
# # stns 2
# lm_fit2 <- lm(formula = ff2, data = data[non_na,] %>% filter(station_id == stns_id[2]))
# summary(lm_fit2)

data <- data %>% mutate(mu = NA , ctemp = NA)
data[non_na & row_id1,] <- data[non_na & row_id1,] %>% 
  mutate(mu = lm_fit1$fitted.values, ctemp = max_temp - mu)
data[non_na & row_id2,] <- data[non_na & row_id2,] %>% 
  mutate(mu = lm_fit2$fitted.values, ctemp = max_temp - mu)

data0 <- data %>% group_by(year, station_name) %>%
  summarise(m_mt = mean(max_temp, na.rm=T),
            m_mu = mean(mu, na.rm=T))
ggplot(data0, aes(x=year, y=m_mt, col=station_name)) +
  geom_line() + geom_point() +
  geom_line(data = data0, aes(y=m_mu))

data0 <- data %>% filter(year %in% 1950:1970) %>%
  group_by(year, month, station_name) %>%
  summarise(m_mt = mean(max_temp, na.rm=T),
            m_mu = mean(mu, na.rm=T))
ggplot(data0, aes(x=year+month/12, y=m_mt, col=station_name)) +
  geom_line() + geom_point() +
  geom_line(data = data0, aes(y=m_mu)) +
  facet_wrap(~station_name)


# nice little plots of the fit: run add infinitum
# some plots might show a break (when the data goes from Ottawa to Toronto)
k <- sample(nrow(data),1)
ii <- k + 1:365
plot(data$time[ii], data$max_temp[ii])
lines(data$time[ii], data$mu[ii], col=2)


# note that this breaks the equalities in the data
# but there is still some heterogeneity in the variance
# N <- 20000
N <- 40000 
plot(data$yday[data$station_id == stns_id[1]][1:N], data$max_temp[1:N], cex=.25)
plot(data$yday[data$station_id == stns_id[1]][1:N], data$ctemp[1:N], cex=.25)

plot(data$yday[data$station_id == stns_id[2]][1:N], data$max_temp[1:N], cex=.25)
plot(data$yday[data$station_id == stns_id[2]][1:N], data$ctemp[1:N], cex=.25)



# Compute pseudo-observations using ECDF ----------------------------------
data <- data %>% mutate(pseudo_temp = NA, diff = NA, w1 = NA, w = NA)
for(id in stns_id){
  cat("Working on ECDF for station ", id, ".\n")
  row_id <- data$station_id == id
  
  for(tt in 1:366){
    if(tt %% 20 == 0) cat("progess: ", round(tt/366*100,2), "%\n")
    
    # construct weights based on yday using gaussian kernel
    data[row_id,] <- data[row_id,] %>%
      mutate(diff = tt - (time %% year_len)) %>%
      mutate(diff = abs(pmin(diff, year_len - diff))) %>%
      mutate(w1 = dnorm(diff,0,1)) %>%
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
plot(data$pseudo_temp[s + 1:100], type="o") # makes sense
plot(data$date[s + 1:1000], data$diff[s + 1:1000]) # linear in time, as expected
plot(data$date[s + 1:1000], data$w[s + 1:1000], type="l") # select only end of year periods

# note that there are equalities, but very very few of them. can disregard.
tab <- table(data[row_id1,]$pseudo_temp)
tab[tab > 1]; sum(tab)
tab <- table(data[row_id2,]$pseudo_temp)
tab[tab > 1]; sum(tab)



# Save for later use ------------------------------------------------------
saveRDS(data, "app/data/data_pp_alt.rds")
