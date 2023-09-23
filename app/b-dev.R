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
stns_id <- c(568)
p <- length(stns_id)
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

data <- data %>% filter(station_id %in% stns_id) %>%
  arrange(station_id, date)

data <- data %>% filter(year %in% 1925:2000)

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
nSC <- 6
s_cols <- paste0("s", 1:nSC)
c_cols <- paste0("c", 1:nSC)
for (i in 1:nSC) data <- data %>%
  mutate(!!s_cols[i] := sin(i*2*pi*time/year_len),
         !!c_cols[i] := cos(i*2*pi*time/year_len))


# create natural spline basis (for long term time trend)
# for this check range for each station
row_ids <- sapply(stns_id, \(id) data$station_id == id)
bss <- lapply(1:p, \(k){
  # data[row_ids[,k],] %>% filter(day == 21, month %in% c(3,6,9,12)) %>% select(time) %>% unlist
  # data[row_ids[,k],] %>% filter(day == 21, month %in% c(6,12)) %>% select(time) %>% unlist
  
  # data[row_ids[,k],] %>% filter(day == 1, month == 1, year %in% seq(1920,2020,20)) %>% select(time) %>% unlist
  data[row_ids[,k],] %>% filter(day == 1, month == 1, year %in% seq(1925,2000,25)) %>% select(time) %>% unlist
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
sc_cols2 <- paste0("(", paste0(c(s_cols, c_cols), collapse = " + "), ")")
# b_cols2 <- paste0("(", paste0(b_cols[[1]], collapse = " + "), ")")
# ff <- paste0("mean_temp ~ ", sc_cols2, "*", b_cols2)
ff <- paste0("mean_temp ~ ", sc_cols2)
lms <- lapply(1:p, \(k){
  # x_cols <- c("mean_temp", s_cols, c_cols, b_cols[[k]])
  # x_cols <- c("mean_temp", s_cols, c_cols)
  # lm(formula = "mean_temp~.", data = data[row_ids[,k] & non_na, x_cols])
  lm(formula = as.formula(ff), data = data[row_ids[,k] & non_na, ])
})

# Bakerville
summary(lms[[1]])

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
data %>% group_by(month, day) %>% summarise(cen = mean(ctemp, na.rm=T), sd = sd(ctemp, na.rm=T))
ggplot(data %>% group_by(month, day) %>% summarise(cen = mean(ctemp, na.rm=T), sd = sd(ctemp, na.rm=T)),
       aes(x=month+day/31, y=cen/sd)) +
  geom_line()

ggplot(data %>% group_by(month, day) %>% reframe(cen = mean(mean_temp, na.rm=T), mu=mean(mu, na.rm=T)),
       aes(x=month+day/31)) +
  geom_line(aes(y=cen)) + geom_line(aes(y=mu), col=2)

data1 <- data %>% filter(year <= 1940) %>% group_by(month, day) %>% reframe(cen = mean(mean_temp, na.rm=T), mu=mean(mu, na.rm=T))
data2 <- data %>% filter(year >= 1940 & year <= 1980) %>% group_by(month, day) %>% reframe(cen = mean(mean_temp, na.rm=T), mu=mean(mu, na.rm=T))
data3 <- data %>% filter(year >= 1980) %>% group_by(month, day) %>% reframe(cen = mean(mean_temp, na.rm=T), mu=mean(mu, na.rm=T))
ggplot(mapping = aes(x=month+day/31)) +
  geom_line(data=data1, aes(y=cen)) + geom_line(data=data1, aes(y=mu)) +
  geom_line(data=data2, aes(y=cen),col=2) + geom_line(data=data2, aes(y=mu), col=2) +
  geom_line(data=data3, aes(y=cen),col=3) + geom_line(data=data3, aes(y=mu), col=3)


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
# saveRDS(data, "app/data/data_pp_debug.rds")
saveRDS(data, "app/data/data_pp_debug.rds")


# plot(diff(Ys[[1]][46355 + (-500:500),1]))

# packages ----------------------------------------------------------------
library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/ms.cpp")
sourceCpp("src/dac_serial.cpp")
source("functions.R")
par(mfrow=c(1,1), mar=c(2,2,1,1))
data <- readRDS("app/data/data_pp_debug.rds")

# load and arrange data ---------------------------------------------------
stns_id <- unique(data$station_id)
stns_name <- unique(data$station_name)

X <- data %>% ungroup %>%
  select(date, station_name, pseudo_temp) %>%
  pivot_wider(names_from = station_name, values_from = pseudo_temp) %>%
  select(any_of(as.character(stns_name))) %>%
  as.matrix()
cities <- colnames(X)

# actual dataset
p <- 2 # number of lags we consider
pp <- p+1 # our datasets will thus have p+1 columns
q <- length(stns_id) # number of stations
Y <- sapply(1:pp, \(r) X[r:(nrow(X)-pp+r),1])
rm(X); gc()
nY <- nrow(Y)

# record which rows to remove when comes the time
na_rows <- which(apply(is.na(Y), 1, any))
n <- nY - length(na_rows)

# function to get specific dataset
pickY <- function(k, YYs = Ys, nna_rows = na_rows){
  Y[-na_rows,] 
}


# Data analysis -----------------------------------------------------------

pw <- T
ss <- seq(2, pp, 1)
if(pw){
  C0 <- matrix(NA, nY, length(ss))
  C0[-na_rows,] <- sapply(ss, \(j) taujack_ms(pickY(0)[,c(1,j)], returnC = T)[,1,2])
}else{
  C0 <- matrix(NA, nY, p)
  C0[-na_rows,] <- dac_serial(pickY(0))
}
p <- length(ss)
pp <- p+1

# pairwise stuff, for comparison
cm0 <- colSums(C0, na.rm=T)/(n*(n-1))
if(pw){
  th0s <- 2*cm0 - 1 # as before...
  g0s <- t(2 * t(C0)/(n-1) - 1 - th0s)
}else{
  th0s <- (2^(1:p)*cm0 - 1)/(2^(1:p) - 1) # as before...
  g0s <- t((2^(1:p)*t(C0)/(n-1) - 1)/(2^(1:p) - 1) - th0s)
}

plot(1:p, th0s, ylim=range(th0s), type="o", cex=.5, pch=19)
abline(h=0, lty=2)
K <- 3000
# zeta0s <- sapply(0:K, \(k) colMeans(g0s[1:(nY-k),]*g0s[(k+1):nY,], na.rm=T)) |> t()
# zeta0s <- sapply(0:K, \(k) sapply(1:ncol(g0s), \(r) cov(g0s[1:(nY-k),r],g0s[(k+1):nY,r], use="pairwise"))) |> t()
zeta0s <- sapply(0:K, \(k) sapply(1:ncol(g0s), \(r) cov(g0s[1:(nY-k),r][-na_rows],g0s[(k+1):nY,r][-na_rows], use="pairwise"))) |> t()
plot(0:K, zeta0s[,1], type="o", xlim=c(0,K))
abline(h=0,lty=2)
for(i in 2:p) lines(0:K, zeta0s[,i], type="o", col=i)
plot(zeta0s[1,1] + 2*cumsum(zeta0s[-1,1]), type="l"); zeta0s[,1]
abline(v=seq(1,365.25*40,365.25), lty=2, col=2)
sig0s <- 4*(zeta0s[1,] + 2*colSums(zeta0s[-1,]))
plot(sig0s, ylim=c(0,max(sig0s)), type="o")
abline(h=0, lty=2)
sqrt(sig0s)
sqrt(n)*th0s
plot(pchisq(n*th0s^2/sig0s, df=1, lower.tail = F))
mean(zeta0s[K - 0:1000,1]); format(mean(zeta0s[K - 0:1000,1]), scientific = F)

       