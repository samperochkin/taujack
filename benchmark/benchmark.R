# File used to perform the benchmarking the algorithms.
# ran with
# R --vanilla --no restore CMD BATCH benchmark/benchmark.R benchmark/log.txt

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
num_rep <- 200
taus <- c(0, .25, .5, .75, 1)
ks <- 7:18 # (n = 2^k)
ps <- seq(2,10,2) # dimensions considered
ps_sub <- ps[c(1,length(ps))] # dimensions considered for bf alg.

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
    
    set.seed(666*s)
    X <- (rnorm(n, 0, sqrt(rho)) + matrix(rnorm(n*10, 0, sqrt(1-rho)), n, 10))
    
    # very fast, needs many replications to get a good estimate
    
    K <- 100
    meds <- summary(microbenchmark(
      "knight_o" = cor.fk(X[,1:2]),
      "knight_e" = taujack_ms(X[,1:2]),
      times = K, unit = "seconds"))$median
    
    time_knight_o <- meds[1]
    time_knight_e <- meds[2]
    
    # still fast, but much less so (unless tau=1)
    time_dac <- as.numeric(rep(NA,length(ps)))
    K <- ifelse(tau == 1 | k <= 12, 50, 5)
    for(r in seq_along(ps)){
      
      if(tau < 1 & k > 16 & ps[r] > 2) next
      
      K <- ifelse(tau == 1 | ps[r] == 2, 50, 10)
      time_dac[r] <- summary(microbenchmark(
        "dac" = taujack_dac(X[,1:ps[r]], thresh=10L),
        times = K, unit = "seconds"))$median
    }
    
    # not fast, unless n is small
    time_bf <- as.numeric(c(NA, NA))
    K <- ifelse(k <= 12, 50, 3)
    if(k <= 16){
      for(r in c(1,2))
        time_bf[r] <- summary(microbenchmark(
          "bf" = taujack_bf(X[,1:ps_sub[r]], seq = T),
          times = K, unit = "seconds"))$median
    }
    
    dt <- data.table(rep_id=rep_id, n=n, tau=tau, rho=rho, p = c(2, 2, ps, ps_sub),
                     fun = c("Knight (original, KO)", "Knight (extended, KE)",
                             rep(c("divide-and-conquer (DAC)"), each=length(ps)),
                             rep(c("brute force (BF)"), each=length(ps_sub))),
                     time = c(time_knight_o, time_knight_e, time_dac, time_bf))
    
    return(dt)
  }, mc.cores = 10) |> rbindlist()
  
  fwrite(times, paste0("benchmark/times", x, ".csv"))
}
