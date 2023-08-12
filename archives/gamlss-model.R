library(tidyverse)


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

table(round(diff(data$time),10))
n <- nrow(data)
data$temp_lag1 <- c(rep(NA,1), data$temp[-n])
data$temp_lag2 <- c(rep(NA,2), data$temp[-((n-1):n)])
data$temp_lag3 <- c(rep(NA,3), data$temp[-((n-2):n)])
data$temp_lag6 <- c(rep(NA,6), data$temp[-((n-5):n)])
data$temp_lag12 <- c(rep(NA,12), data$temp[-((n-11):n)])
data$temp_lag24 <- c(rep(NA,24), data$temp[-((n-23):n)])
data <- data %>% filter(!is.na(temp),
                        !is.na(data$temp_lag1),
                        !is.na(data$temp_lag2),
                        !is.na(data$temp_lag3),
                        !is.na(data$temp_lag6),
                        !is.na(data$temp_lag12),
                        !is.na(data$temp_lag24))
data <- data[order(data$station_id),]

year_len <- 365 + 6/24 + 9/60/24 + 9/60^2/24
moon_len <- 1 + 0/24 + 50/60/24

data <- data %>% mutate(s1_time = sin(2*pi*time/year_len),
                        c1_time = cos(2*pi*time/year_len),
                        # s2_time = sin(4*pi*time/year_len),
                        # c2_time = cos(4*pi*time/year_len),
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


# data0 <- data[1:25000, c("temp", "time", "hour",
data0 <- data[, c("temp", "date", "time", "month", "day", "hour",
                         grep("temp_lag", names(data), value = T),
                          grep("_time", names(data), value = T),
                         grep("_moon", names(data), value = T),
                         grep("_hour", names(data), value = T))]
library(gamlss)
# ff <- temp ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12 + temp_lag24)*
ff <- temp ~ (temp_lag1*temp_lag2)*
  (s1_time + c1_time)*
  # (s1_time + c1_time + s2_time + c2_time)*
  (s1_hour + c1_hour + s2_hour + c2_hour +
     s3_hour + c3_hour + s4_hour + c4_hour +
     s5_hour + c5_hour + s6_hour + c6_hour +
     s7_hour + c7_hour)
mod <- gamlss(formula = ff, ff, data = data0)



# data0 <- data[data$year <= 1980 & data$station_id == stns_id[1],]
# mod <- readRDS("mods.rds")[[as.character(stns_id[1])]]

mod
plot(mod$residuals[1:5000])
hist(mod$residuals, breaks=50)


plot(mod$mu.fv[1:1000], type="o", cex=.25)
plot(mod$sigma.fv[1:1000], type="o", cex=.25)
plot(mod$mu.fv[1:1000]-data0$temp[1:1000], type="o", cex=.25)

# data0 <- data[data$year <= 1980 & data$station_id == stns_id[1],]
res <- (data0$temp - mod$mu.fv)/mod$sigma.fv
qqnorm(res)
res <- (data0$temp - mod$mu.fv)/mod$sigma.fv


ind <- which(data0$month == 7 &  data0$day %in% 1:7 & data0$hour == 12)
ind <- which(data0$month == 7 &  data0$day %in% 20:27 & data0$hour == 12)
ind <- which(data0$month == 2 &  data0$day %in% 1:7 & data0$hour == 7)
data00 <- data0[ind,]

hist((data00$temp - mod$mu.fv[ind])/mod$sigma.fv[ind])
plot(mod$mu.fv[1:1000], type="o", cex=.25)

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
abline(v=seq(0,24*10000,24), col=2, lty=2)
acf(x, 72, ylim=c(-.05,.05), na.action = na.pass)
abline(v=seq(0,24*10000,24), col=2, lty=2)

pacf(x, 500*3, ylim=c(-.05,.05), na.action = na.pass)
acf(x, 500*3, ylim=c(-.05,.05), na.action = na.pass)

pacf(x, 500*10, ylim=c(-.05,.05), na.action = na.pass)
acf(x, 500*10, ylim=c(-.05,.05), na.action = na.pass)



