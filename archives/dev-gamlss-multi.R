library(Rcpp)
library(tidyverse)
sourceCpp("src/dac.cpp")
sourceCpp("src/ms.cpp")


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
data <- data %>% mutate(time = as.integer(data$date) + data$hour/24,
                        yday = lubridate::yday(date),
                        week = lubridate::week(date),
                        seas = lunar::terrestrial.season(date))
str(data)

getLag <- function(dat){
  print(table(round(diff(dat$time),10)))
  n <- nrow(dat)
  dat$temp_lag1 <- c(rep(NA,1), dat$temp[-n])
  dat$temp_lag2 <- c(rep(NA,2), dat$temp[-((n-1):n)])
  dat$temp_lag3 <- c(rep(NA,3), dat$temp[-((n-2):n)])
  dat$temp_lag6 <- c(rep(NA,6), dat$temp[-((n-5):n)])
  dat$temp_lag12 <- c(rep(NA,12), dat$temp[-((n-11):n)])
  dat$temp_lag24 <- c(rep(NA,24), dat$temp[-((n-23):n)])
  dat <- dat %>% filter(!is.na(temp),
                          !is.na(dat$temp_lag1),
                          !is.na(dat$temp_lag2),
                          !is.na(dat$temp_lag3),
                          !is.na(dat$temp_lag6),
                          !is.na(dat$temp_lag12),
                          !is.na(dat$temp_lag24))
  return(dat)
}


data <- lapply(unique(data$station_id),
               \(st) getLag(data[data$station_id == st,])) |>
    do.call(what = "rbind")



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


data0 <- data[, c("station_id", "temp", "time", "hour", "seas",
                         grep("temp_lag", names(data), value = T),
                         grep("_time", names(data), value = T),
                         grep("_moon", names(data), value = T),
                         grep("_hour", names(data), value = T))]

data0 <- data[data$year <= 1980, c("station_id", "date", "temp", "time", "hour", "seas",
                  grep("temp_lag", names(data), value = T),
                  grep("_time", names(data), value = T),
                  grep("_moon", names(data), value = T),
                  grep("_hour", names(data), value = T))]


library(gamlss)
ff <- temp ~ (temp_lag1*temp_lag2 + temp_lag6 + temp_lag12 + temp_lag24)*
  (s1_time + c1_time + s2_time + c2_time)*
  (s1_moon + c1_moon + s1_hour + c1_hour + s2_hour + c2_hour)

# data0$res <- 0
# mods <- list()
# for(st in unique(data0$station_id)){
#   mod <- gamlss(formula = ff, ff, data = data0[data0$station_id == st,])
#   mods[[as.character(st)]] <- mod
#   data0[data0$station_id == st,]$res <- 
#     (data0[data0$station_id == st,]$temp - mod$mu.fv)/mod$sigma.fv
# }
# saveRDS(mods, "mods.rds")

mods <- readRDS("mods.rds")
data0$res <- 0
for(st in unique(data0$station_id)){
  mod <- mods[[as.character(st)]]
  data0[data0$station_id == st,]$res <-
    (data0[data0$station_id == st,]$temp - mod$mu.fv)/mod$sigma.fv
}



###
names(data)
data00 <- data0 %>% dplyr::select(date, seas, hour, station_id, res) %>%
  pivot_wider(names_from = station_id,
              values_from = res) %>% 
  drop_na()

# X <- data00 %>% ungroup %>% 
X <- data00 %>% ungroup %>% filter(seas == "Winter") %>%
  dplyr::select(any_of(as.character(stns_id))) %>%
  as.matrix()

pairs(X[1:5000,], cex=.1)
U <- apply(X,2,copula::pobs)
pairs(U[1:5000,], cex=.1)
Th <- pcaPP::cor.fk(U)
Th
nrow(X)






# season comparison -------------------------------------------------------
Xs <- lapply(unique(data00$seas), \(ss){
  data00 %>% ungroup %>% filter(seas == ss) %>%
    dplyr::select(any_of(as.character(stns_id))) %>%
    as.matrix()
})
names(Xs) <- unique(data00$seas)


# pairwise tau 
Ths <- lapply(Xs, \(X){
  n <- nrow(X)
  d <- ncol(X)
  M <- array(n-1, dim=c(d,d,n))
  for(i in 1:(d-1)){
    X <- X[order(X[,i]),]
    for(j in (i+1):d){
      M[i,j,] <- M[j,i,] <- n - 1 - countSwaps(X[,j])
    }
  }
  
  Th <- apply(2*M/(n-1)-1, c(1,2), mean)
  g <- 2*apply(combn(d,2), 2, \(ij) M[ij[1],ij[2],])/(n-1)-1
  S <- 4*var(g)/n
    
  list("Tau" = Th, "tau" = Th[t(combn(d,2))], "S" = S)
})
lapply(Ths, "[[", "Tau")


# multivariate concordance for NS provinces.
Ths <- lapply(Xs, dac)
sapply(Ths, \(C){
  n <- length(C)
  CC <- (2^3 * C/(n-1) - 1)/(2^3 - 1)
  c(mean(CC),
    4*var(CC)/n)
})

sapply(Ths, \(C){
  n <- length(C)
  CC <- (2^3 * C/(n-1) - 1)/(2^3 - 1)
  c(mean(CC),
    4*var(CC)/n)
})



