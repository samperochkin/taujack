library(Rcpp)
sourceCpp("src/ms.cpp")
sourceCpp("src/bf.cpp")
sourceCpp("src/bfs.cpp")
sourceCpp("src/dac_serial.cpp")

# quick sanity check with Knight's extended alg.
d <- 2
n <- 1000
rho = sin(.5 * pi/2)
X <- rnorm(n*d,0,sqrt(1-rho)) |> matrix(nrow=n, ncol=d) + rnorm(n,0,sqrt(rho))
C0 <- n-1-countSwaps(X[order(X[,1]),2])[rank(X[,1])]
C1 <- bruteForceGen(X, serial=T)
C2 <- dac_serial(X, 25)
all(C0 == C1)
all(C0 == C2)


# quick sanity check with (basic) bruteForce
d <- 6
n <- 1000
rho = sin(.5 * pi/2)
X <- rnorm(n*d,0,sqrt(1-rho)) |> matrix(nrow=n, ncol=d) + rnorm(n,0,sqrt(rho))
C0 <- sapply(2:d, \(k) bruteForce(X[,1:k]))
C1 <- bruteForceGen(X, serial=T)
C2 <- dac_serial(X, 25)
all(C0 == C1)
all(C0 == C2)

