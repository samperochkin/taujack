library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/dac_serial.cpp")
# sourceCpp("src/dac_serial_debug.cpp")


d <- 6
m <- 5
n <- m^d

N <- 100
runtimes <- replicate(N, {
  # co-monotone
  X <- replicate(d, 1:n)
  t1 <- system.time(dac_serial(X, 25))[3]
  
  # perfectly random
  X <- replicate(d, 1:m, simplify = F) |> expand.grid() |> as.matrix() + runif(n*d,-.01,.01)
  # X <- X[sample(n),]
  t2 <- system.time(dac_serial(X, 25))[3]

  # random
  X <- replicate(d, sample(n))
  t3 <- system.time(dac_serial(X, 25))[3]
  
  # dependent (rho = 1/2)
  X <- rnorm(n*d) |> matrix(nrow=n, ncol=d) + rnorm(n)
  t4 <- system.time(dac_serial(X, 25))[3]
  
  # dependent (tau = .5)
  rho = sin(.5 * pi/2)
  X <- rnorm(n*d,0,sqrt(1-rho)) |> matrix(nrow=n, ncol=d) + rnorm(n,0,sqrt(rho))
  t5 <- system.time(dac_serial(X, 25))[3]
  
  # X shaped
  W <- sample(0:1,n*d,TRUE) |> matrix(nrow=n, ncol=d)
  X <- W*replicate(d, 1:n) + (1-W)*replicate(d, n:1) + rnorm(n*d,0,n/10)
  # pairs(X[sample(n,1000),], pch=19, cex=.25)
  # U <- apply(X,2,rank)
  # pairs(U[sample(n,1000),], pch=19, cex=.25)
  t6 <- system.time(dac_serial(X, 25))[3]
  
  c("co-monotone" = t1,
    "perf. uniform" = t2, "uniform" = t3,
    "normal (rho=.5)" = t4, "normal (tau=.5)" = t5,
    "X shaped" = t6)
})

runtimes
rowMeans(runtimes)
apply(runtimes,1,quantile)





d <- 2
n <- 7000

rho = sin(.5 * pi/2)
X <- rnorm(n*d,0,sqrt(1-rho)) |> matrix(nrow=n, ncol=d) + rnorm(n,0,sqrt(rho))

system.time(dac_serial(X, 25))[3]
system.time(bruteForce(X))[3]
system.time(cor(X, method="kendall"))[3]
