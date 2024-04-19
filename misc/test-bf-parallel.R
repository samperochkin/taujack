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

K <- 10
# X <- mvtnorm::rmvnorm(5e3, rep(0,5), S[1:5,1:5])

rho <- sin(pi*.5/2)
X <- mvtnorm::rmvnorm(2^18, rep(0,K), (1-rho)*diag(K) + rho)

system.time(taujack_dac(X[,1:K],10))
system.time(taujack_bf(X[,1:K],T))
system.time(taujack_bf_parallel(X[,1:5],T,mc_cores = 3))
system.time(taujack_dac_parallel(X[,1:5],T,mc_cores = 3))

system.time(taujack_dac(X[,1:4]))
system.time({
  taujack_dac(X[,1:4])
  taujack_dac(X[,2:5])
  taujack_dac(X[,c(3:5,1)])
  taujack_dac(X[,c(4:5,1:2)])
  taujack_dac(X[,c(5,1:3)])
})
#   taujack_dac(X[,c(1,3,4)])
#   taujack_dac(X[,c(1,4,5)])
# })

