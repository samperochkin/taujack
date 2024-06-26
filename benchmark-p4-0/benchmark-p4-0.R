# File used to perform the benchmarking the algorithms.
# ran with
# R --vanilla --no restore CMD BATCH benchmark-p4-0/benchmark-p4-0.R benchmark-p4-0/log.txt

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
num_rep <- 1
taus <- 0
ks <- 7:25 # (n = 2^k)
p <- c(4) # dimensions considered

sim_grid <- expand.grid(rep_id = 1:num_rep, tau = taus, k = ks)
sim_grid$n <- 2^sim_grid$k
sim_grid$rho <- sin(sim_grid$tau*pi/2)
cat("Number of row in sim_grid: ", nrow(sim_grid), "\n")

ns <- 2^ks
K <- 5

mclapply(1:nrow(sim_grid), \(s){
  
  cat("------------------------- s=", s, "\n")

    cat(".\n")
    
    rep_id <- sim_grid[s,]$rep_id
    n <- sim_grid[s,]$n
    k <- sim_grid[s,]$k
    tau <- sim_grid[s,]$tau
    rho <- sim_grid[s,]$rho
    
    set.seed(34238*s)
    if(rho >= 0){
      X <- (rnorm(n, 0, sqrt(rho)) + matrix(rnorm(n*p, 0, sqrt(1-rho)), n, p))
    }else{
      X <- mvtnorm::rmvnorm(n, sigma = (1-rho)*diag(p) + rho)
    }
    
    # still fast, but much less so (unless tau=1)
    time_dac <- summary(microbenchmark(
      "dac" = taujack_dac(X[,1:p], thresh=as.integer(2^9)),
      times = K, unit = "seconds"))$median

    dt <- data.table(rep_id=rep_id, n=n, tau=tau, rho=rho, p = p,
                     fun = "divide-and-conquer (DAC)",
                     time = time_dac)
    
    fwrite(dt, paste0("benchmark-p4-0/times", s, ".csv"))
    NULL
}, mc.cores = 6) |> rbindlist()
