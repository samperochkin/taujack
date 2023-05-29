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
taus <- c(0, .5, 1)
ks <- 8:18 # (n = 2^k)
ps <- seq(2,10,2) # dimensions considered

sim_grid <- expand.grid(rep_id = 1:num_rep, tau = taus, k = ks)
sim_grid$n <- 2^sim_grid$k
sim_grid$rho <- sin(sim_grid$tau*pi/2)

ns <- 2^ks

for(r in seq_along(ns)){
  times <- mclapply(which(sim_grid$n == ns[r]), \(s){
    
    cat(".\n")
    
    rep_id <- sim_grid[s,]$rep_id
    n <- sim_grid[s,]$n
    tau <- sim_grid[s,]$tau
    rho <- sim_grid[s,]$rho
  
    set.seed(34*s)
    X <- (rnorm(n, 0, sqrt(rho)) + matrix(rnorm(n*15, 0, sqrt(1-rho)), n, 15))
  
    # very fast, needs many replications to get a good estimate
    time_knight_o <- system.time(replicate(100, cor.fk(X[,1:2])))[[3]]/100
    time_knight_e <- system.time(replicate(100, taujack_ms(X[,1:2])))[[3]]/100
    
    # fast, but not too fast (unless tau=1)
    time_dac <- numeric(length(ps))
    for(r in seq_along(ps)){
      K <- ifelse(tau == 1, 50, 5)
      time_dac[r] <- system.time(replicate(K, taujack_dac(X[,1:ps[r]], thresh=100L)))[[3]]/K
    }
  
    # not fast, unless n is small
    time_bf <- numeric(2)
    ps_sub <- ps[c(1,length(ps))]
    for(r in c(1,2)){
      K <- ifelse(n < 2^10, 10, 1)
      time_bf[r] <- system.time(replicate(K, taujack_bf(X[,1:ps_sub[r]])))[[3]]/K
    }
    
    dt <- data.table(rep_id=rep_id, n=n, tau=tau, rho=rho, p = c(2, 2, ps, ps_sub),
                     fun = c("Knight (original)", "Knight (extended)",
                             rep(c("d-a-c"), each=length(ps)),
                             rep(c("brute force"), each=length(ps_sub))),
                     time = c(time_knight_o, time_knight_e, time_dac, time_bf))
    
  
    return(dt)
  }, mc.cores = 12) |> rbindlist()




  fwrite(times, paste0("times", r, ".csv"))
}