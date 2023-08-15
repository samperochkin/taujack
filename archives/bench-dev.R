# File used to perform the benchmarking the algorithms.

# Packages ----------------------------------------------------------------
library(Rcpp)
library(data.table)
library(ggplot2)
library(parallel)

# functions to benchmark
library(pcaPP) # cor.fk
sourceCpp("src/ms.cpp")   # for Knight's extended alg.
sourceCpp("src/dac.cpp")  # divide-and-conquer alg.
sourceCpp("src/bf.cpp")   # brute force alg.
source("functions.R")       # wrappers (performs re-ordering if necessary)


# Benchmark ---------------------------------------------------------------
num_rep <- 100
taus <- c(0)
ks <- 8:18 # (n = 2^k)
# ps <- seq(2,10,2) # dimensions considered

sim_grid <- expand.grid(k = ks, tau = taus, rep_id = 1:num_rep)
sim_grid$n <- 2^sim_grid$k
sim_grid$rho <- sin(sim_grid$tau*pi/2)

times <- mclapply(sample(nrow(sim_grid)), \(s){
  
  cat(".\n")
  
  rep_id <- sim_grid[s,]$rep_id
  n <- sim_grid[s,]$n
  tau <- sim_grid[s,]$tau
  rho <- sim_grid[s,]$rho
  
  set.seed((66+99)*s)
  X <- (rnorm(n, 0, sqrt(rho)) + matrix(rnorm(n*15, 0, sqrt(1-rho)), n, 15))
  
  # fast, but not too fast (unless tau=1)
  K <- 5
  time_dac <- system.time(replicate(K, taujack_dac(X[,1:8], thresh=100L)))[[3]]/K

  dt <- data.table(rep_id=rep_id, n=n, tau=tau, rho=rho, p = c(8),
                   fun = "d-a-c", time = time_dac)
  
  return(dt)
}, mc.cores = 12) |> rbindlist()




fwrite(times, "times4.csv")
