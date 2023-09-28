

###########################################################################
# Main script (data analysis of daily mean temperatures) ------------------
###########################################################################

# packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(Rcpp)
sourceCpp("src/bf.cpp") # naive/brute force
sourceCpp("src/ms.cpp") # knight's extended
sourceCpp("src/dac_seq.cpp") # divide-and-conquer
source("functions.R") # R wrapper functions and more



# load and arrange data ---------------------------------------------------
data <- readRDS("app/data/data_pp.rds")
stns_id <- unique(data$station_id)
stns_name <- unique(data$station_name)

# reorder to have Toronto in middle
stns_id <- stns_id[match(c("Welland", "Toronto", "Ottawa"), stns_name)]
stns_name <- stns_name[match(c("Welland", "Toronto", "Ottawa"), stns_name)]

# check range of data
data %>% group_by(station_name) %>% 
  filter(!is.na(pseudo_temp)) %>%
  summarize(min_date = min(date), max_date = max(date), n = length(date))

# data <- data %>% mutate(pseudo_temp = ifelse(month %in% c(12,1,2), pseudo_temp, NA))

# create dataset of univariate time series
X1 <- data %>% ungroup %>%
  select(date, month, station_name, pseudo_temp) %>%
  filter(month %in% 5:10) %>%
  pivot_wider(names_from = station_name, values_from = pseudo_temp) %>%
  select(any_of(as.character(stns_name))) %>%
  as.matrix()

X2 <- data %>% ungroup %>%
  select(date, month, station_name, pseudo_temp) %>%
  filter(month == c(12,1,2,3)) %>%
  pivot_wider(names_from = station_name, values_from = pseudo_temp) %>%
  select(any_of(as.character(stns_name))) %>%
  as.matrix()


N <- 50000
L <- 10 # number of lags we consider
p <- L+1 # our datasets will thus have p+1 columns
q <- 1 # number of stations

z <- rnorm(N + L)
z[1:N] <- z[1:N + L] + .25 * z[1:N + L - 1] + .1 * z[1:N + L - 2] + .01 * z[1:N + L - 3] + .1 * rnorm(N)
X1 <- sapply(0:L, \(k) z[1:N + L - k])
X2 <- sapply(0:L, \(k) z[1:N + k] + rnorm(N))
Ys <- list(X1, X2)
nY <- nrow(Ys[[1]])

# record which rows to remove when comes the time
na_rowss <- lapply(Ys, \(Y) which(apply(is.na(Y[,1:p]), 1, any)))
ns <- nY - sapply(na_rowss, length)

# function to get specific dataset (relies on global env)
pickY <- function(k) Ys[[k]][-na_rowss[[k]],] 

# round(pcaPP::cor.fk(pickY(1)),3)
# round(pcaPP::cor.fk(pickY(2)),3)
round(pcaPP::cor.fk(Ys[[1]]),3)
round(pcaPP::cor.fk(Ys[[2]]),3)

# pairwise stuff, for comparison
K <- 500
res_pw <- lapply(1:q, \(k){
  cat("Work on k=", k,".\n")
  # C <- sapply(2:p, \(r) dac_seq(pickY(k)[,c(1,r)])) |>
  #   C2Cna(na_rows = na_rowss[[k]]) |> Cna2tsz(K_serial = p, pw = T)
  C <- sapply(2:p, \(r) dac_seq(Ys[[k]][,c(1,r)])) |> Cna2tsz(K_serial = K, pw = T)
})

# check zetas
k <- 1
Zh <- res_pw[[k]]$Zh
Zh <- Zh + aperm(Zh, c(2,1,3))
Zhcs <- 2*apply(Zh, c(1,2), \(z) z[1] + 2*cumsum(z[-1])) |> aperm(perm = c(2,3,1))
i <- 1
par(mfrow = c(2,5), mar=c(2,2,3,1))
for(j in 1:min(L,10)){plot(Zhcs[i,j,], main = paste0("pair: (",i,",",j,")"), type="l", ylim=c(-.3,.75)); abline(h=0, lty=2, col=2)}

Ths <- sapply(res_pw, "[[", "Th"); Ths
par(mfrow = c(1,2), mar=c(2,2,3,1))
plot(Ths[,1], type="l"); abline(h=0, lty=2, col=2)
plot(Ths[,2], type="l"); abline(h=0, lty=2, col=2)
sapply(lapply(res_pw, "[[", "Sh"), matrixcalc::is.positive.definite)

# check zetas
k <- 1
Zh <- res_pw[[k]]$Zh
Zh <- Zh + aperm(Zh, c(2,1,3))
Zhcs <- 2*apply(Zh, c(1,2), \(z) z[1] + 2*cumsum(z[-1])) |> aperm(perm = c(2,3,1))

# zcs <- 4*(Zh[1,1,1] + cumsum((N-1:K)/N * Zh[1,1,-1]))
# plot(zcs, type="l")
# zcs <- 4*(Zh[1,1,1] + cumsum(Zh[1,1,-1]))
# plot(zcs, type="l")
# res <- replicate(50, {
#   L <- 10
#   p <- L+1
#   print("a")
#   z <- rnorm(N + L)
#   z[1:N] <- z[1:N + L] + .25 * z[1:N + L - 1] + .1 * z[1:N + L - 2] + .01 * z[1:N + L - 3] + .1 * rnorm(N)
#   X1 <- sapply(0:L, \(k) z[1:N + L - k])
#   (sapply(2:p, \(r) dac_seq(X1[,c(1,r)])) |> Cna2tsz(K_serial = 0, pw = T))$Th
# })
# N*cov(t(res))


i <- 1
par(mfrow = c(2,5), mar=c(2,2,3,1))
for(j in 1:min(L,10)){plot(Zh[i,j,], main = paste0("pair: (",i,",",j,")"), type="l"); abline(h=0, lty=2, col=2)}
mean(Zh[1,4,-(1:20)])
for(j in 1:min(L,10)){plot(Zhcs[i,j,], main = paste0("pair: (",i,",",j,")"), type="l", ylim=c(-.1,.5)); abline(h=0, lty=2, col=2)}
Shs <- lapply(lapply(res_pw, "[[", "Zh"), Zh2Sh, K=15)
sapply(Shs, \(S) matrixcalc::is.positive.definite(S[1:15,1:15]))

