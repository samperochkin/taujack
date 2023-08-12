
# packages ----------------------------------------------------------------
library(gamlss)
library(tidyverse)



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
                        c3_time = cos(6*pi*time/year_len),
                        s1_hour = sin(2*pi*time),
                        c1_hour = cos(2*pi*time),
                        s2_hour = sin(4*pi*time),
                        c2_hour = cos(4*pi*time),
                        s3_hour = sin(6*pi*time),
                        c3_hour = cos(6*pi*time),
                        s4_hour = sin(8*pi*time),
                        c4_hour = cos(8*pi*time),
                        s5_hour = sin(10*pi*time),
                        c5_hour = cos(10*pi*time),
                        s6_hour = sin(12*pi*time),
                        c6_hour = cos(12*pi*time),
                        s7_hour = sin(14*pi*time),
                        c7_hour = cos(14*pi*time))

# subset data so that no NA creeps up in usused columns
data0 <- data[, c("keep", "temp", "date", "time", "month", "day", "hour",
                  grep("temp_lag", names(data), value = T),
                  grep("_time", names(data), value = T),
                  grep("_hour", names(data), value = T))]
ind <- which(data0$keep)

ff <- temp ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12 + temp_lag24)*
  (s1_time + c1_time + s2_time + c2_time + s3_time + c3_time)

mclapply(0:23, \(h){
  # nams <- paste0("fv", c("m", "s", "n"), h)
  indh <- ind[data0[ind,]$hour == h]
  
  mod <- gamlss(formula = ff, ff, ff,
                family = gamlss.dist::SN1(),
                data = data0[indh, 1:33],
                control = gamlss.control(c.crit = .01, n.cyc = 150),
                i.control = glim.control(cc = .01))
  
  saveRDS(cbind(mod$mu.coefficients,
                mod$sigma.coefficients,
                mod$nu.coefficients), paste0("app/P_", h, ".rds"))
})

# 
# 
# for(h in 0:23){
#   cat("\n")
#   cat("----------------------------------\n")
#   cat("Fitting model four hour == ", h, "\n")
#   
#   nams <- paste0("fv", c("m", "s", "n"), h)
#   indh <- ind[data0[ind,]$hour == h]
#   
#   if(h == 0){
#     mod <- gamlss(formula = ff, ff, ff,
#                   family = gamlss.dist::SN1(),
#                   data = data0[indh, 1:33],
#                   control = gamlss.control(c.crit = .01, n.cyc = 150),
#                   i.control = glim.control(cc = .01))
#     Pm <- mod$mu.coefficients
#     Ps <- mod$sigma.coefficients
#     Pn <- mod$nu.coefficients
#   }else{
#     mod <- gamlss(formula = ff, ff, ff,
#                   family = gamlss.dist::SN2(),
#                   mu.start = mean(mod$mu.fv),
#                   sigma.start = mean(mod$sigma.fv),
#                   nu.start = mean(mod$nu.fv),
#                   data = data0[indh, 1:33],
#                   control = gamlss.control(c.crit = .01, n.cyc = 150),
#                   i.control = glim.control(cc = .01))
#     Pm <- rbind(Pm, mod$mu.coefficients)
#     Ps <- rbind(Ps, mod$sigma.coefficients)
#     Pn <- rbind(Pn, mod$nu.coefficients)
#   }
#   
#   for(nam in nams) data0[[nam]] <- as.numeric(NA)
#   data0[[nams[1]]][indh] <- mod$mu.fv
#   data0[[nams[2]]][indh] <- mod$sigma.fv
#   data0[[nams[3]]][indh] <- mod$nu.fv
# }
# 




ff1 <- temp - fv2 - fv3 - fv4 - fv5 - fv6 - fv7 - fv8 - fv9 - fv10 ~ temp_lag1*temp_lag2 + temp_lag6 + temp_lag12 +
  s1_time + c1_time + s2_time + c2_time + s3_time + c3_time +
  s1_hour + c1_hour + s2_hour + c2_hour +
  s3_hour + c3_hour + s4_hour + c4_hour +
  s5_hour + c5_hour + s6_hour + c6_hour +
  s7_hour + c7_hour

ff2 <- temp - fv1 - fv3 - fv4 - fv5 - fv6 - fv7 - fv8 - fv9 - fv10 ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12):
  (s1_time + c1_time + s2_time + c2_time + s3_time + c3_time)

ff3 <- temp - fv1 - fv2 - fv4 - fv5 - fv6 - fv7 - fv8 - fv9 - fv10 ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12):
  (s1_hour + c1_hour + s2_hour + c2_hour + s3_hour + c3_hour + s4_hour + c4_hour)

ff4 <- temp - fv1 - fv2 - fv3 - fv5 - fv6 - fv7 - fv8 - fv9 - fv10 ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12):
  (s5_hour + c5_hour + s6_hour + c6_hour + s7_hour + c7_hour)

ff5 <- temp - fv1 - fv2 - fv3 - fv4 - fv6 - fv7 - fv8 - fv9 - fv10 ~ (s1_time + c1_time + s2_time + c2_time + s3_time + c3_time):
  (s1_hour + c1_hour + s2_hour + c2_hour + s3_hour + c3_hour + s4_hour + c4_hour)

ff6 <- temp - fv1 - fv2 - fv3 - fv4 - fv5 - fv7 - fv8 - fv9 - fv10 ~ (s1_time + c1_time + s2_time + c2_time + s3_time + c3_time):
  (s5_hour + c5_hour + s6_hour + c6_hour + s7_hour + c7_hour)

ff7 <- temp - fv1 - fv2 - fv3 - fv4 - fv5 - fv6 - fv8 - fv9 - fv10 ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12):
  (s1_time + c1_time + s2_time + c2_time + s3_time + c3_time):
  (s1_hour + c1_hour + s2_hour + c2_hour)

ff8 <- temp - fv1 - fv2 - fv3 - fv4 - fv5 - fv6 - fv7 - fv9 - fv10 ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12):
  (s1_time + c1_time + s2_time + c2_time + s3_time + c3_time):
  (s3_hour + c3_hour + s4_hour + c4_hour)

ff9 <- temp - fv1 - fv2 - fv3 - fv4 - fv5 - fv6 - fv7 - fv8 - fv10 ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12):
  (s1_time + c1_time + s2_time + c2_time + s3_time + c3_time):
  (s5_hour + c5_hour + s6_hour)

ff10 <- temp - fv1 - fv2 - fv3 - fv4 - fv5 - fv6 - fv7 - fv8 - fv9 ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12):
  (s1_time + c1_time + s2_time + c2_time + s3_time + c3_time):
  (c6_hour + s7_hour + c7_hour)

# ff0 <- temp ~ temp_lag1*temp_lag2 + temp_lag6 + temp_lag12 +
#   s1_time + c1_time + s2_time + c2_time + s3_time + c3_time +
#   s1_hour + c1_hour + s2_hour + c2_hour +
#   s3_hour + c3_hour + s4_hour + c4_hour +
#   s5_hour + c5_hour + s6_hour + c6_hour +
#   s7_hour + c7_hour
ff0 <- temp ~ temp_lag1*temp_lag2 + temp_lag6 + temp_lag12 +
  s1_time + c1_time + s2_time + c2_time + s3_time + c3_time

ffs <- list(ff1,ff2,ff3,ff4,ff5,ff6,ff7,ff8,ff9,ff10)

data0$fv1 <- as.numeric(0)
data0$fv2 <- as.numeric(0)
data0$fv3 <- as.numeric(0)
data0$fv4 <- as.numeric(0)
data0$fv5 <- as.numeric(0)
data0$fv6 <- as.numeric(0)
data0$fv7 <- as.numeric(0)
data0$fv8 <- as.numeric(0)
data0$fv9 <- as.numeric(0)
data0$fv10 <- as.numeric(0)



# Initial run with subset of data -----------------------------------------

D_diff <- Inf
counter <- 0
while(D_diff > 1){
  counter <- counter + 1
  cat("\n",
      "--------------------\n",
      "--------------------\n",
      "iteration number ", counter, "\n")
  
  if(counter == 1){
    D_old <- Inf
  }else{
    D_old <- D_new
  }
  for(i in 1:8){
    cat("\n", "--------------------\n", "model ", i, "\n")
    mod <- gamlss(formula = ffs[[i]], ff0, ff0, family = gamlss.dist::SN2(), data = data0[ind,],
                  mu.start = data0[ind,paste0("fvm",i)][[1]],
                  sigma.start = data0[ind,paste0("fvs",i)][[1]],
                  nu.start = data0[ind,paste0("fvn",i)][[1]],
                  control = gamlss.control(c.crit = .1), i.control = glim.control(cc = .1))
    data0[ind,paste0("fvm",i)][[1]] <- mod$mu.fv
    data0[ind,paste0("fvs",i)][[1]] <- mod$sigma.fv
    data0[ind,paste0("fvn",i)][[1]] <- mod$nu.fv
  }
  D_new <- mod$G.deviance
  D_diff <- D_old - D_new
  cat("\n", "--------------------\n", "D_diff: ", D_diff, "\n")
}





D_diff <- Inf
counter <- 0
while(D_diff > 1){
  counter <- counter + 1
  cat("\n",
      "--------------------\n",
      "--------------------\n",
      "iteration number ", counter, "\n")

  if(counter == 1){
    D_old <- Inf
  }else{
    D_old <- D_new
  }
  for(i in 1:8){
    cat("\n", "--------------------\n", "model ", i, "\n")
    mod <- gamlss(formula = ffs[[i]], ff0, ff0, family = gamlss.dist::SN2(), data = data0[ind,],
                  mu.start = data0[ind,paste0("fvm",i)][[1]],
                  sigma.start = data0[ind,paste0("fvs",i)][[1]],
                  nu.start = data0[ind,paste0("fvn",i)][[1]],
                  control = gamlss.control(c.crit = .1), i.control = glim.control(cc = .1))
    data0[ind,paste0("fvm",i)][[1]] <- mod$mu.fv
    data0[ind,paste0("fvs",i)][[1]] <- mod$sigma.fv
    data0[ind,paste0("fvn",i)][[1]] <- mod$nu.fv
  }
  D_new <- mod$G.deviance
  D_diff <- D_old - D_new
  cat("\n", "--------------------\n", "D_diff: ", D_diff, "\n")
}

saveRDS(data0, "data_w_fv.rds")

data0$fvm <- as.numeric(0)
data0$fvs <- as.numeric(0)
data0$fvn <- as.numeric(0)
fvm <- rowSums(sapply(1:10, \(i) data0[ind,paste0("fvm",i)][[1]]))
fvs <- rowSums(sapply(1:10, \(i) data0[ind,paste0("fvs",i)][[1]]))
fvn <- rowSums(sapply(1:10, \(i) data0[ind,paste0("fvn",i)][[1]]))
data0$res <- as.numeric(0)
data0[ind,]$res <- pSN2(data0[ind,]$temp, fvm, fvs, fvn)
hist(data0[ind,]$res, breaks=100)

for(m in 12:1){
  data0$keep2 <- as.logical(NA)
  data0$keep2 <- data0$keep & (data0$day %in% 1:7 & data0$month == m)
  ind2 <- which(data0$keep2)
  hist(data0[ind2,]$res, breaks=40, xlim = c(-10,10), ylim = c(0,1), probability = T)
}


qqnorm(sort(data0$res[1:10000]), type="l")
abline(0,1, col=2, lty=2)


# let's compute smoothed empirical cdf for the normalized residuals
# and further transform them
new_res <- rep(NA, length(res))
my_grid <- expand.grid(month=1:12, day=1:31, hour=0:23)
for(i in 1:nrow(my_grid)){
  if(i %% 100 == 0) print(i/nrow(my_grid))
  
  m <- grid[i,"month"]
  d <- grid[i,"day"]
  h <- grid[i,"hour"]
  
  ind <- which(data0$day == d & data0$month == m)
  if(length(ind) == 0) next
  
  ind <- which(floor(data0$time) %in% unlist(lapply(ind, \(i) seq(data0$date[i] - 7, data0$date[i] + 7,1))) &
                 data0$hour == h)
  
  # hist(res[ind], breaks=50, probability = T)
  ran <- c(floor(min(res[ind])),ceiling(max(res[ind])))
  e <- round(0.3 * diff(ran))
  den <- density(res[ind], bw=.2, from=ran[1]-e, to=ran[2]+e)
  # lines(den$x, den$y)
  
  cdf <- cumsum(den$y*(den$x[2]-den$x[1]))
  cdf <- cdf/cdf[length(cdf)]
  
  dd <- data.frame(x = den$x, y = cdf)
  loess_mod <- loess(y~x, data=dd)
  
  new_res[ind] <- predict(loess_mod, res[ind])
}

saveRDS(new_res, "new_res.rds")


par(mfrow=c(1,1))
xid <- data0$time*24
xid <- xid - min(xid) + 1
x <- numeric(diff(range(xid)))
# xid <- xid[data0$seas == "Summer"]
# x[xid] <- data0$temp_res[data0$seas == "Summer"]
x[xid] <- new_res
x[-xid] <- NA
which(is.na(x))

n <- length(x)
pacf(x, 72, ylim=c(-.05,.05), na.action = na.pass)
abline(v=seq(0,24*1000,6), col=4, lty=2)
abline(v=seq(0,24*1000,12), col=3, lty=2)
abline(v=seq(0,24*1000,24), col=2, lty=2)
acf(x, 72, ylim=c(-.05,.05), na.action = na.pass)
abline(v=seq(0,24*1000,6), col=4, lty=2)
abline(v=seq(0,24*1000,12), col=3, lty=2)
abline(v=seq(0,24*1000,24), col=2, lty=2)

pacf(x, 500*3, ylim=c(-.05,.05), na.action = na.pass)
acf(x, 500*3, ylim=c(-.05,.05), na.action = na.pass)
acf(x, 500*3, ylim=c(-.05,.05), xlim = c(500,600), na.action = na.pass)

pacf(x, 500*10, ylim=c(-.05,.05), na.action = na.pass)
acf(x, 500*10, ylim=c(-.05,.05), na.action = na.pass)



plot(x[1:5000], x[100 + 1:5000])
xx <- rank(x)/(length(x)+1)
plot(xx[1:5000], xx[100 + 1:5000])


