# File used to perform the benchmarking the algorithms.
# ran with
# R --vanilla --no restore CMD BATCH benchmark/benchmark.R benchmark/log.txt

# Packages ----------------------------------------------------------------
library(Rcpp)
library(data.table)
library(ggplot2)
library(parallel)

# functions to benchmark
library(pcaPP) # cor.fk
# sourceCpp("src/ms.cpp")   # for Knight's extended alg.
sourceCpp("src/dac_seq.cpp")  # divide-and-conquer alg.
# sourceCpp("src/bf.cpp")   # brute force alg.
source("functions.R")       # wrappers (performs re-ordering if necessary)


# Benchmark ---------------------------------------------------------------
num_rep <- 100
taus <- c(0, .25, .5, .75, 1)
ks <- 8:19 # (n = 2^k)
ps <- seq(2,10,2) # dimensions considered
ps_sub <- ps[c(1,length(ps))] # dimensions considered for bf alg.

sim_grid <- expand.grid(rep_id = 1:num_rep, tau = taus, k = ks)
sim_grid$n <- 2^sim_grid$k
sim_grid$rho <- sin(sim_grid$tau*pi/2)
cat("Number of row in sim_grid: ", nrow(sim_grid), "\n")

ns <- 2^ks
x <- 6
s <- which(sim_grid$k == ks[x])[1]
# times <- lapply(which(sim_grid$k == ks[x])[-(1:400)], \(s){
  
  cat(s, "\n")
  
  rep_id <- sim_grid[s,]$rep_id
  n <- sim_grid[s,]$n
  k <- sim_grid[s,]$k
  tau <- sim_grid[s,]$tau
  rho <- sim_grid[s,]$rho
  
  set.seed(34*s)
  X <- rnorm(n, 0, sqrt(rho)) + matrix(rnorm(n*10, 0, sqrt(1-rho)), n, 10)
  
  # very fast, needs many replications to get a good estimate
  time_knight_o <- system.time(replicate(20, cor.fk(X[,1:2])))[[3]]/20
  time_knight_e <- system.time(replicate(10, taujack_ms(X[,1:2])))[[3]]/10
  
  # still fast, but much less so (unless tau=1)
  time_dac <- as.numeric(rep(NA,length(ps)))
  K <- ifelse(tau == 1 | k <= 12, 50, 5)
  r <- 2
  # if(!(tau < 1 & k > 16)){
  #   for(r in seq_along(ps))
        taujack_dac(X[,1:ps[r]], thresh=25L)
        time_dac[r] <- system.time(replicate(K, taujack_dac(X[,1:ps[r]], thresh=25L)))[[3]]/K
  # }
  
  # not fast, unless n is small
  # time_bf <- as.numeric(c(NA, NA))
  # K <- ifelse(k <= 12, 10, 1)
  # if(k <= 16){
  #   for(r in c(1,2))
  #     time_bf[r] <- system.time(replicate(K, taujack_bf(X[,1:ps_sub[r]], seq = T)))[[3]]/K
  # }
  # 
  # dt <- data.table(rep_id=rep_id, n=n, tau=tau, rho=rho, p = c(2, 2, ps, ps_sub),
  #                  fun = c("Knight (original, KO)", "Knight (extended, KE)",
  #                          rep(c("divide-and-conquer (DAC)"), each=length(ps)),
  #                          rep(c("brute force (BF)"), each=length(ps_sub))),
  #                  time = c(time_knight_o, time_knight_e, time_dac, time_bf))
  
