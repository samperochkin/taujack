library(parallel)
library(Rcpp)
sourceCpp("../src/bf.cpp")   # brute force alg.

print("Started from the bottom. Now we here.")

seq <- T
n <- 50
X <- replicate(5, sample(n))
mc_cores <- 5

nn <- round(seq(1,n,length.out=mc_cores+1))
nn <- matrix(c(nn[1], sapply(nn[-c(1,mc_cores+1)], \(x) c(x,x+1)), nn[mc_cores+1]),
             mc_cores, 2, byrow = T)

clus <- parallel::makeCluster(mc_cores)
parallel::clusterExport(clus, list("seq"), envir = environment())
parallel::clusterEvalQ(clus, {
  library(taujackmininal)
})

parallel::parLapply(cl = clus,
                    X = apply(nn, 1, \(x) X[seq(x[1],x[2]),], simplify = F),
                    fun = \(X_sub) bruteForceR(X_sub))

stopCluster(clus)
