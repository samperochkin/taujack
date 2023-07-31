library(tidyverse)


filepaths <- list.files("data", pattern = "temp_data", full.names = T)
stns_id <- readRDS("data/stns_id.rds")

filepaths <- filepaths[grepl(stns_id[1], filepaths) |
                         grepl(stns_id[2], filepaths) |
                         grepl(stns_id[3], filepaths) |
                         grepl(stns_id[4], filepaths)]
data <- lapply(filepaths, readRDS) |> dplyr::bind_rows()
data <- data %>% mutate(date = as.Date(date),
                        year = as.integer(year),
                        month = as.integer(month),
                        day = as.integer(day),
                        hour = lubridate::hour(strptime(hour, format = '%H:%M')))
data <- data %>% mutate(time = as.integer(data$date) + data$hour/24)
str(data)

data <- data %>% filter(!is.na(temp))
data <- data[order(data$station_id),]

my_pobs = function(temp, bw=5) {
  # Extend range of density estimate beyond data
  ran <- c(floor(min(temp)),ceiling(max(temp)))
  e <- round(0.3 * diff(ran))
  
  # dens <- ecdf(temp)
  dens <- density(temp, bw=bw, from=ran[1]-e, to=ran[2]+e)
  
  # Fit a model
  # mod <- loess(cdf~x, data=dd)
  splinefun(x=cumsum(dens$y), method="hyman")(temp)
}

year_len <- 365 + 6/24 + 9/60/24 + 9/60^2/24
# moon_len <- 27 + 7/24 + 43/60/24
# moon_len <- 29
moon_len <- 1 + 0/24 + 50/60/24

rr <- function(time, date, hour, temp, ky = 4, kd = 3, cross = T, kl = 2){
  # rr <- function(time, temp){
  B <- cbind(1,
             sapply(1:ky, \(k) sin(2*k*pi*time/year_len)),
             sapply(1:ky, \(k) cos(2*k*pi*time/year_len)),
             sapply(1:kd, \(k) sin(2*k*pi*time)),
             sapply(1:kd, \(k) cos(2*k*pi*time)))
  
  if(cross){
    l1 <- ((2*ky)+1)
    l2 <- (2*kd)
    B <- c(list(B), lapply(2:l1, \(b) b*B[,(l1+1):(l1+l2)])) |> do.call(what="cbind")
  }

  # phase <- lunar::lunar.phase(date, hour)
  if(kl > 0){
    B <- cbind(B,
               # sapply(1:kl, \(k) sin(k*phase)),
               # sapply(1:kl, \(k) cos(k*phase)))
               sapply(1:kl, \(k) sin(2*k*time/moon_len)),
               sapply(1:kl, \(k) cos(2*k*time/moon_len)))
  }

  m <- B %*% (MASS::ginv(B) %*% temp)
  # my_pobs(temp-m)
  temp-m
}

data <- data %>% 
  group_by(station_id) %>%
  mutate(temp_res = rr(time, date, hour, temp, ky = 3, kd = 3, cross=F, kl = 1)) %>% 
  group_by(station_id, month, day, hour) %>%
  mutate(temp_res = my_pobs(temp_res))


plot(data$temp_res[1:5000], type="o", cex=.25)

par(mfrow=c(1,1))
data0 <- data %>% filter(station_id == stns_id[1])
n <- nrow(data0)
pacf(data0$temp_res[seq(1,n,24)], 500*30, ylim=c(-.05,.05))
acf(data0$temp_res[seq(1,n,24)], 500*30, ylim=c(-.05,.05))
aa <- acf(data0$temp_res[seq(1,n,24)], 500*30, ylim=c(-.05,.05))
plot(data0$hour[seq(1,n,24)], c(aa$acf), ylim=c(-.05,.05))

acf(data0$temp_res[seq(1,n,24)], 2000, ylim=c(-.05,.05), col=1:24)
acf(data0$temp_res[seq(1,n,24)], 2000, xlim=c(700,700+24*32), ylim=c(-.05,.05), col=1:24)



names(data)
data0 <- data %>% dplyr::select(date, time, year, month, day, hour, station_id, temp_res) %>%
  pivot_wider(names_from = station_id,
              values_from = temp_res) %>% 
  drop_na()

X <- data0 %>% ungroup %>% 
  dplyr::select(any_of(as.character(stns_id))) %>%
  as.matrix()


data0 <- data0 %>% mutate(seas = lunar::terrestrial.season(date))
lapply(unique(data0$seas), \(s){
  pcaPP::cor.fk(X[data0$seas == s,])
})
lapply(unique(data0$seas), \(s){
  pairs(apply(X[data0$seas == s,][1:1000,],2, copula::pobs), cex=.1)
})

data0 <- data0 %>% group_by(hour) %>% mutate(seas_nd = paste0(seas, "_", hour < 12))
data0$seas_nd[1:10]

lapply(unique(data0[data0$seas == "Summer",]$seas_nd), \(s){
  print(s)
  print(pcaPP::cor.fk(X[data0$seas_nd == s,]))
  NULL
})

data0 <- data0 %>% group_by(hour) %>% mutate(seas_h = paste0(seas, "_", hour))
data0$seas_nd[1:10]

lapply(unique(data0[data0$seas == "Summer",]$seas_h), \(s){
  print(s)
  print(pcaPP::cor.fk(X[data0$seas_h == s,]))
  NULL
})

Th <- sapply(unique(data0[data0$seas == "Summer",]$seas_h), \(s){
  pcaPP::cor.fk(X[data0$seas_h == s,])
}, simplify = "array")

par(mfrow = c(5,6), mar=c(2,1,2,2))
for(i in 1:7){
  for(j in (i+1):8){
    plot(Th[i,j,], main = paste0(i, "_", j), type="o", cex=.25, ylim=c(0,1))
  }
}
