# File used to perform the benchmarking the algorithms.
# ran with
# R --vanilla --no restore CMD BATCH benchmark-p4-unif/benchmark-p4-unif.R benchmark-p4-unif/log.txt

# Packages ----------------------------------------------------------------
library(Rcpp)
library(data.table)
library(ggplot2)
library(parallel)
library(microbenchmark, lib.loc = "~/R/x86_64-pc-linux-gnu-library/4.3/")

# functions to benchmark
library(pcaPP) # cor.fk
sourceCpp("src/dac_seq.cpp")  # divide-and-conquer alg.
sourceCpp("src/dac_seq_tune.cpp")  # divide-and-conquer alg.
source("functions.R")       # wrappers (performs re-ordering if necessary)


# Benchmark ---------------------------------------------------------------
num_rep <- 1
taus <- 0
n0s <- 5:34 # (n = 2^k)
p <- c(4) # dimensions considered

sim_grid <- expand.grid(rep_id = 1:num_rep, n0 = n0s, p=p)
sim_grid$n <- sim_grid$n0^p
cat("Number of row in sim_grid: ", nrow(sim_grid), "\n")

K <- 25
mclapply(1:nrow(sim_grid), \(s){
    
  cat("------------------------- n=", n0s[x], "\n")

  rep_id <- sim_grid[s,]$rep_id
  n0 <- sim_grid[s,]$n0

  X <- expand.grid(1:n0,1:n0,1:n0,1:n0) |> as.matrix() #+ (1:n0^4)/(n0^4+1)
  n <- nrow(X)

  gc()
  time_dac <- summary(microbenchmark(
    "dac" = taujack_dac_tune(X, thresh=as.integer(2^9), tune=3L),
    times = K, unit = "seconds"))$median

  dt <- data.table(rep_id=rep_id, n=n, p = p,
                   fun = "divide-and-conquer (DAC)",
                   time = time_dac)
    
  fwrite(dt, paste0("benchmark-p4-unif/times", s, ".csv"))
}, mc.cores = 10)
