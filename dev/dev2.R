
# ii <- which.max(abs(g[,1]))
# plot(Ys[[1]][ii + (-500:500),1])
# plot(diff(Ys[[1]][ii + (-500:500),1]))

rr <- 5000
plot((ii + (-rr:rr)) / 365,data[data$station_id == stns_id[1],]$mean_temp[ii + (-rr:rr)], pch=19, cex=.25)
lines((ii + (-rr:rr)) / 365,data[data$station_id == stns_id[1],]$mu[ii + (-rr:rr)], col=2, lwd=2)
abline(v=1:200, lty=2)
plot((ii + (-rr:rr)) / 365,data[data$station_id == stns_id[1],]$ctemp[ii + (-rr:rr)], pch=19, cex=.25)
abline(v=1:200, lty=2)
plot((ii + (-rr:rr)) / 365,data[data$station_id == stns_id[1],]$pseudo_temp[ii + (-rr:rr)], pch=19, cex=.25)
abline(v=1:200, lty=2)
plot((ii + (-rr:(rr-1))) / 365,diff(data[data$station_id == stns_id[1],]$pseudo_temp[ii + (-rr:rr)]), pch=19, cex=.25)
abline(v=1:200, lty=2)

plot(g[1:(nrow(g)-100),1],g[101:nrow(g),1])
plot(g[1:(nrow(g)-1),1],g[2:nrow(g),1])
hist(g[,4])


# packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/ms.cpp")
sourceCpp("src/dac_serial.cpp")
source("functions.R")

# load and arrange data ---------------------------------------------------
data <- readRDS("app/data/data_pp.rds")
# data <- data %>% mutate(seas = lunar::terrestrial.season(date))
stns_id <- unique(data$station_id)
stns_name <- unique(data$station_name)

X <- data %>% ungroup %>%
  select(date, station_name, pseudo_temp) %>%
  pivot_wider(names_from = station_name, values_from = pseudo_temp) %>%
  select(any_of(as.character(stns_name))) %>%
  as.matrix()
cities <- colnames(X)

X <- replicate(q, rnorm(nrow(X)))
X[5:nrow(X),] <- X[5:nrow(X),] + .5 * X[4:(nrow(X)-1),] + .25 * X[1:(nrow(X)-4),] + rnorm(nrow(X)-4,0,.01)

# actual dataset
p <- 5 # number of lags we consider
pp <- p+1 # our datasets will thus have p+1 columns
q <- length(stns_id) # number of stations
Ys <- lapply(1:q, \(k) sapply(1:pp, \(r) X[r:(nrow(X)-pp+r),k]))
rm(X); gc()
nY <- nrow(Ys[[1]])
ns <- rep(nY, q)


Ys[[1]] <- mvtnorm::rmvnorm(nY, rep(0,pp), .75*diag(pp) + .25)


# Data analysis -----------------------------------------------------------

############################################
#### Multivariate tau for each stations ####
############################################

# concordance matrices - keep dims of Y to ensure equally spaced timestamps
C <- dac_serial(Ys[[1]])
cm <- colSums(C, na.rm=T)/(ns[1]*(ns[1]-1))
tt <- (2^(1:p) * cm - 1)/(2^(1:p) - 1)
g <- t((2^(1:p) * t(C)/(ns[1]-1) - 1)/(2^(1:p) - 1) - tt)
plot(1:p, tt, ylim=c(0,max(tt)), type="o", col=1, lwd=1.5)

# compute the K first terms of the sum in Sen's Theorem 1
K <- 200
zetas <- sapply(0:K, \(k) colMeans(g[1:(nY-k),]*g[(k+1):nY,], na.rm=T)) |> t()
# plot(0:K, zetas[,1], type="o", xlim=c(0,K))
# for(i in 2:p) lines(0:K, zetas[,i], type="o", col=i)


# K <- 1000
# r <- 10
# # zz <- colMeans(g[1:(nY-K),]*g[(K+1):nY,], na.rm=T)
# zz <- sapply(0:K, \(k) mean(g[1:(nY-k),r]*g[(k+1):nY,r], na.rm=T))
# plot(0:K/365.25, zz)
# abline(v=1:10, lty=2)
# plot(0:K/365.25, cumsum(zz))
# abline(v=1:10, lty=2)

# corresponding jackknife variance
sigs <- 4*(zetas[1,] + 2*colSums(zetas[-1,]))
# plot(sigs, ylim=c(0,max(sigs)), type="o")
lines(sigs, ylim=c(0,max(sigs)), type="o")
abline(h=0, lty=2)










######################

p <- 10 # number of lags we consider
pp <- p+1 # our datasets will thus have p+1 columns
N <- nY/6
samp <- replicate(100, {
  X <- rnorm(N)
  X[5:N] <- X[5:N] + .5 * X[4:(N-1)] + .25 * X[1:(N-4)] + rnorm(N-4,0,.01)
  Y <- sapply(1:pp, \(r) X[r:(N-pp+r)])
  C <- dac_serial(Y)
  cm <- colSums(C, na.rm=T)/(N*(N-1))
  tt <- (2^(1:p) * cm - 1)/(2^(1:p) - 1)
})

samp <- t(samp - rowMeans(samp))
pairs(samp[,1:6])
