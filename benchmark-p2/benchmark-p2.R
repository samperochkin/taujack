# File used to perform the benchmarking the algorithms.
# ran with
# R --vanilla --no restore CMD BATCH benchmark-p2/benchmark-p2.R benchmark-p2/log.txt

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
num_rep <- 5
taus <- 0
ks <- 25 # (n = 2^k)
p <- c(2) # dimensions considered

sim_grid <- expand.grid(rep_id = 1:num_rep, tau = taus, k = ks)
sim_grid$n <- 2^sim_grid$k
sim_grid$rho <- sin(sim_grid$tau*pi/2)
cat("Number of row in sim_grid: ", nrow(sim_grid), "\n")

ns <- 2^ks

for(x in seq_along(ks)){
  cat("------------------------- k=", ks[x], "\n")
  times <- mclapply(which(sim_grid$k == ks[x]), \(s){
    
    cat(".\n")
    
    rep_id <- sim_grid[s,]$rep_id
    n <- sim_grid[s,]$n
    k <- sim_grid[s,]$k
    tau <- sim_grid[s,]$tau
    rho <- sim_grid[s,]$rho
    
    set.seed(53424*s + 1)
    X <- (rnorm(n, 0, sqrt(rho)) + matrix(rnorm(n*p, 0, sqrt(1-rho)), n, p))
    K <- 1
    
    # still fast, but much less so (unless tau=1)
    time_dac <- summary(microbenchmark(
      "dac" = taujack_dac(X[,1:p], thresh=as.integer(2^9)),
      times = K, unit = "seconds"))$median

    dt <- data.table(rep_id=rep_id, n=n, tau=tau, rho=rho, p = p,
                     fun = "divide-and-conquer (DAC)",
                     time = time_dac)
    
    return(dt)
  }, mc.cores = 15) |> rbindlist()
  
  fwrite(times, paste0("benchmark-p2/times_0_29_", 20+x, ".csv"))
}
