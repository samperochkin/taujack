# Packages ----------------------------------------------------------------
library(Rcpp)
library(data.table)
library(ggplot2)
library(parallel)
library(microbenchmark, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.3/")

# functions to benchmark
library(pcaPP) # cor.fk
sourceCpp("src/ms.cpp")   # for Knight's extended alg.
sourceCpp("src/dac_seq.cpp")  # divide-and-conquer alg.
sourceCpp("src/bf.cpp")   # brute force alg.
source("functions.R")       # wrappers (performs re-ordering if necessary)


# Benchmark ---------------------------------------------------------------
n <- 3e4
k <- 25
ps <- c(2,6)
rho <- .5

S <- diag(k)
S <- S + matrix(runif(k^2,-.2,.6),k,k)
S <- (S + t(S))/2

S <- diag(k)
a <- .9
for(r in 1:(k-1)) S[cbind(1:(k-r),(r+1):k)] <- a^r
S <- (S + t(S)) - diag(k)
image(S)
eigen(S)$val

# X <- (rnorm(n, 0, sqrt(rho)) + matrix(rnorm(n*k, 0, sqrt(1-rho)), n, k))
X <- mvtnorm::rmvnorm(n, rep(0,k), S)
Th <- cor.fk(X)
hist(Th)

K <- 10
summary(microbenchmark(
  "dac" = taujack_dac(X, thresh=10L),
  "bf" = taujack_bf(X, seq=T),
  times = K, unit = "seconds"))$median

ij <- which(Th == min(Th), arr.ind = T)[,1]
X0 <- X[,c(ij, (1:k)[-ij])]
summary(microbenchmark(
  "dac" = taujack_dac(X0, thresh=10L),
  "bf" = taujack_bf(X0, seq=T),
  times = K, unit = "seconds"))$median


ij <- order(rowSums(Th))
X0 <- X[,ij]
summary(microbenchmark(
  "dac" = taujack_dac(X0, thresh=10L),
  "bf" = taujack_bf(X0, seq=T),
  times = K, unit = "seconds"))$median



n <- 6e4
X <- mvtnorm::rmvnorm(n, rep(0,k), S)
Th <- cor.fk(X)

ij <- which(Th == min(Th), arr.ind = T)[1,]
for(r in 4:k){
  ij <- c(ij,(1:k)[-ij][which(Th[ij,-ij] == min(Th[ij,-ij]), arr.ind = T)[1,2]])
}
ij <- c(ij, (1:k)[-ij])
X0 <- X[,ij]
summary(microbenchmark(
  "dac" = taujack_dac(X0, thresh=10L),
  "bf" = taujack_bf(X0, seq=T),
  times = K, unit = "seconds"))$median

X <- mvtnorm::rmvnorm(6e4, rep(0,k), S)

K <- 3
res <- sapply(seq(5e3,10*5e3, 5e3), \(n){
  print(n)
  Th <- cor.fk(X[1:n,])
  
  ij <- which(Th == min(Th), arr.ind = T)[1,]
  for(r in 4:k){
    ij <- c(ij,(1:k)[-ij][which(Th[ij,-ij] == min(Th[ij,-ij]), arr.ind = T)[1,2]])
  }
  ij <- c(ij, (1:k)[-ij])
  X0 <- X[1:n,ij]
  summary(microbenchmark(
    "dac" = taujack_dac(X0, thresh=10L),
    times = K, unit = "seconds"))$median
})

plot(res)

q <- 25
plot(log(seq(1,2e24, length.out=2000),2)^q, type="l")
lines(seq(1,2e24, length.out=2000)^2, col=2)

for(q in seq(50,100,10)){
  plot(log(seq(1,2^q, length.out=2000),2)^q, type="l")
  lines(seq(1,2^q, length.out=2000)^2, col=2)
}
