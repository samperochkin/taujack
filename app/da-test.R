
# plot(diff(Ys[[1]][46355 + (-500:500),1]))

# packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/ms.cpp")
sourceCpp("src/dac_serial.cpp")
source("functions.R")
par(mfrow=c(1,1), mar=c(2,2,1,1))

# load and arrange data ---------------------------------------------------
data <- readRDS("app/data/data_pp.rds")
# data <- readRDS("app/data/data_pp1.rds")
# data <- readRDS("app/data/data_pp2.rds")
# data <- data %>% mutate(seas = lunar::terrestrial.season(date))
stns_id <- unique(data$station_id)
stns_name <- unique(data$station_name)

X0 <- data %>% ungroup %>%
  select(date, station_name, pseudo_temp) %>%
  pivot_wider(names_from = station_name, values_from = pseudo_temp) %>%
  select(any_of(as.character(stns_name))) %>%
  as.matrix()
cities <- colnames(X0)

# actual dataset
p <- 10 # number of lags we consider
d <- pp <- p+1 # our datasets will thus have p+1 columns
X0 <- sapply(1:pp, \(r) X0[r:(nrow(X0)-pp+r),1])
# ids <- which(!apply(is.na(X),1,any))
# n <- length(ids)
# X[!is.na(X)] <- X[!is.na(X)] + rnorm(sum(!is.na(X)),0,1e-5)

R <- 7500
ii <- which(!sapply(seq(100,nrow(X0)-R,100), \(x) any(is.na(X0[x + 1:R,]))))[5]*100
# ii <- 7600 + 7500
ii <- 20000
R <- 5000
X <- X0[ii + 1:R,]
plot(is.na(X[,1]), type="o"); sum(is.na(X[,1]))
# X <- do.call("rbind", replicate(100, X[ii + 1:R,], simplify = F))
# X <- X + rnorm(n*pp, 0 ,1)
# n <- nrow(X)
# ids <- 1:nrow(X)
ids <- which(apply(!is.na(X), 1, all))
n <- length(ids)

C <- matrix(NA, nrow(X), ncol(X)-1)
C[ids,] <- dac_serial(X[ids,])
# C <- bruteForce(X)

hist(C, breaks=20)

cs <- colSums(C, na.rm=T)
C_jack <- t((2^(2:d -1) * 2*(cs/2 - t(C))/(n*(n-1)) - 1)/(2^(2:d -1) - 1))
hist(C_jack[,1], breaks=100)
tt <- (2^(2:d - 1) * cs/(n*(n-1)) - 1)/(2^(2:d - 1) - 1)
tt
plot(tt)

g <- t((2^(2:d - 1) * t(C)/(n-1) - 1)/(2^(2:d - 1) - 1) - tt)
K <- 2000
zetas <- as.matrix(sapply(0:K, \(k) colMeans(g[1:(n-k),,drop=F]*g[(k+1):n,,drop=F], na.rm=T)))
if(d > 2) zetas <- t(zetas)
plot(zetas[,1]); abline(h=0, lty=2, col=2)
plot(zetas[1,1] + cumsum(zetas[-1,1]))
abline(h=0, lty=2)
plot(diff(zetas[1,1] + 2*cumsum(zetas[-1,1])))
plot(rev(cummean(rev(zetas[-1,1]))))
sigs <- 4*(zetas[1,] + 2*colSums(zetas[-1,,drop=F]))
sigs



K <- 1000
ks <- seq(10,K,10)
res <- sapply(ks, \(i){
  zetas <- as.matrix(sapply(0:i, \(k) colMeans(g[1:(n-k),,drop=F]*g[(k+1):n,,drop=F], na.rm=T)))
  if(d > 2) zetas <- t(zetas)
  sigs <- 4*(zetas[1,] + 2*colSums(zetas[-1,,drop=F]))
  print(round(sigs,2))
  sigs
})

plot(ks, res[9,])
