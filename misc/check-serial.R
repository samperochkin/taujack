# Script used in developing the version for serial dependence 


# Packages ----------------------------------------------------------------
library(Rcpp)
sourceCpp("src/dac_serial.cpp")
sourceCpp("src/bf.cpp")


# Validation of results ---------------------------------------------------
n <- 500
d <- 5
X <- matrix(runif(n*d),n,d)

C0 <- sapply(2:d, \(k) bruteForce(X[,1:k]))
C1 <- dac_serial(X, thresh = 10)
all(C0 - C1 == 0)



# Compare runtimes --------------------------------------------------------
n <- 30000
# n <- 50000
d <- 6
X <- matrix(runif(n*d),n,d)
system.time(bruteForce(X))
system.time(dac_serial(X, thresh = 25))



# Validate estimators -----------------------------------------------------

# function to introduce serial dependence
sDep <- function(X, K){
  n <- nrow(X)
  X[-(1:K),] <- t(sapply((K+1):n, \(i){
    colSums(X[i - (0:K),])
  }))
  X
}

# setup
n <- 1000
N <- n*600
d <- 4

# generate data with classical dependence
X0 <- matrix(rnorm(N*d),N,d) + rnorm(n,0,.25)

# introduce serial dependence
K <- 5
X0 <- sDep(X0, K)

# scale between 0 and 1 and plot
X0 <- apply(X0,2,rank)/(N+1)
plot(X[1:min(n,100),1], type="o")

# approximate "true" values (th0 and sig0)
C0 <- dac_serial(X0, thresh = 25, brute_force = F)[,d-1]
th0 <- (2^(d-1) * sum(C0)/(N*(N-1)) - 1)/(2^(d-1) - 1)
g0 <- (2^(d-1)*C0/(N-1) - 1)/(2^(d-1) - 1) - th0
ss0 <- sapply(0:15, \(k){
  mean(g0[1:(N-k)]*g0[(k+1):N])
})
plot(0:15, ss0, main="Summand in the variance formula of Sen")
abline(h=0, lty=2, col=2)
abline(v=K+1, lty=2, col=3)
# sig0 <- 2*(ss0[1] + 2*sum(ss0[-1]))
sig0 <- ss0[1] + 2*sum(ss0[2:(K+1)])

# compute estimators for N/n (very nearly) independent datasets
ths <- NULL
sigs <- NULL
for(k in 1:(N/n)){
  X <- X0[(k-1)*n + 1:n,]
  C <- dac_serial(X, thresh = 100, brute_force = F)[,d-1]
  th <- (2^(d-1) * sum(C)/(n*(n-1)) - 1)/(2^(d-1) - 1)
  ths <- c(ths, th)
  g <- (2^(d-1)*C/(n-1) - 1)/(2^(d-1) - 1) - th
  
  ss <- sapply(0:15, \(k){
    mean(g[1:(n-k)]*g[(k+1):n])
  })
  # plot(ss)
  # abline(h=0, lty=2, col=2)
  
  # s <- 2*(ss[1] + 2*sum(ss[-1]))
  s <- ss[1] + 2*sum(ss[2:(K+1)])
  sigs <- c(sigs, s)
}

# histogram to check if sig0 matches what we observed
hist(sqrt(n)*(ths-th0), probability = T, breaks=13)
lines(seq(-10,10,.01),dnorm(seq(-10,10,.01),0,2*sqrt(sig0)))

# histogram of sigma values just for fun
hist(sigs)
abline(v=sig0, lty=2, col=2)

# three quantites that should match
4*sig0
4*mean(sigs)
n*var(ths)

