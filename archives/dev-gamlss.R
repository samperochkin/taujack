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
                        s2_time = sin(4*pi*time/year_len),
                        c2_time = cos(4*pi*time/year_len),
                        s1_moon = sin(2*pi*time/moon_len),
                        c1_moon = cos(2*pi*time/moon_len),
                        s1_hour = sin(2*pi*time),
                        c1_hour = cos(2*pi*time),
                        s2_hour = sin(4*pi*time),
                        c2_hour = cos(4*pi*time))


# data0 <- data[1:25000, c("temp", "time", "hour",
data0 <- data[, c("temp", "time", "hour",
                         grep("temp_lag", names(data), value = T),
                          grep("_time", names(data), value = T),
                         grep("_moon", names(data), value = T),
                         grep("_hour", names(data), value = T))]
library(gamlss)
# ff <- temp ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12 + temp_lag24)*
#   (s1_time + c1_time + s2_time + c2_time)*
#   (s1_moon + c1_moon + s1_hour + c1_hour + s2_hour + c2_hour)
# mod <- gamlss(formula = ff, ff, data = data0)

data0 <- data[data$year <= 1980 & data$station_id == stns_id[1],]

mod <- readRDS("mods.rds")[[as.character(stns_id[1])]]

mod
plot(mod$residuals)
hist(mod$residuals, breaks=50)


plot(mod$mu.fv[1:1000], type="o", cex=.25)
plot(mod$sigma.fv[1:1000], type="o", cex=.25)
plot(mod$mu.fv[1:1000]-data0$temp[1:1000], type="o", cex=.25)

# data0 <- data[data$year <= 1980 & data$station_id == stns_id[1],]
res <- (data0$temp - mod$mu.fv)/mod$sigma.fv
qqnorm(res)

a


xx <- res[1:50000]
n <- length(xx)
names(xx) <- NULL
spe <- spectrum(xx, log="no", plot=F)
spx <- spe$freq/2*pi*24*365.25
spy <- 2*spe$spec
plot(spy~spx,xlab="cycles/year",ylab="spectral density",type="l", xaxt="n")
axis(1, at = spx[seq(1,n/2,length.out=5)],
     labels = round(1/spx[seq(1,n/2,length.out=5)]))


plot(spy~spx,xlab="cycles/year",ylab="spectral density",type="l", xaxt="n", xlim = c(0,5))
axis(1, at = spx[seq(1,n/2,length.out=5)], labels = round(1/spx[seq(1,n/2,length.out=5)]))
spx[which.max(spy)]

spe

n <- length(spe$freq)
plot(spe$freq, spe$spec, cex=.25)

plot(1/spe$freq[1:24], spe$spec[1:8760], cex=.25)

plot(log(1/spe$freq[1:8760]/8760), spe$spec[1:8760], cex=.25)
which.max(spe$spec[1:8760])/8760


par(mfrow=c(1,1))
xid <- data0$time*24
xid <- xid - min(xid) + 1
x <- numeric(diff(range(xid)))
# xid <- xid[data0$seas == "Summer"]
# x[xid] <- data0$temp_res[data0$seas == "Summer"]
x[xid] <- res
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
acf(x, 500*30, ylim=c(-.05,.05), na.action = na.pass)



