library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/bf_full.cpp")
sourceCpp("src/bfs.cpp")
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

c(dac_serial(X, 25)) == bruteForce(X)

system.time(dac_serial(X, 25))[3]
system.time(bruteForce(X))[3]
system.time(cor(X, method="kendall"))[3]






d <- 10
n <- 10000

rho = sin(.5 * pi/2)
X <- rnorm(n*d,0,sqrt(1-rho)) |> matrix(nrow=n, ncol=d) + rnorm(n,0,sqrt(rho))

bruteForceFull(X)

system.time(bruteForce(X))[3]
system.time(bruteForceGen(X, serial = T))[3]
system.time(bruteForceGen(X, all = T))[3]


profmem::profmem(bruteForceGen(X))
profmem::profmem(bruteForceGen(X, serial = T))
profmem::profmem(bruteForceGen(X, all = T))



d <- 60
n <- 30000

k <- 4
x <- rnorm(n+d+k)
y <- numeric(n+d)
b <- seq(1,0,length.out=k+1)[-(k+1)]^2
for(i in 1:n) y[i] <- sum(b*x[k + i - k:1])
y <- y + rnorm(n+d,0,.1)
plot(y[1:100], type="o")
X <- sapply(1:d, \(r) y[r:(n+r)])

res1 <- sapply(2:d, \(k){
  print(k)
  system.time(dac_serial(X[,1:k], 25))[[3]]
})
res2 <- sapply(2:d, \(k){
  print(k)
  system.time(bruteForceGen(X[,1:k], serial = T))[[3]]
})
plot(res1, type="o", ylim = range(c(res1,res2)))
lines(res2, type="o")

mean(diff(res1)[(d-10):d - 2])
mean(diff(res2)[(d-10):d - 2])

system.time(dac_serial(X, 25))[3]
system.time(bruteForce(X))[3]
system.time(bruteForceGen(X, serial = T))[3]









d <- 2
q <- 20
n <- 30000

k <- 4
x1 <- rnorm(n+q+k)
x2 <- rnorm(n+q+k) + x1

y1 <- numeric(n+q)
y2 <- numeric(n+q)

b1 <- seq(1,0,length.out=k+1)[-(k+1)]^1.5
b2 <- seq(1,0,length.out=k+1)[-(k+1)]^2

for(i in 1:n){
  y1[i] <- sum(b1*x1[k + i - k:1])
  y2[i] <- sum(b2*x2[k + i - k:1])
} 
z <- rnorm(n+q,0,.1)
y1 <- y1 + z
y2 <- y2 + rnorm(n+q,0,.1) + .5*z
plot(y1[1:100], type="o")
lines(y2[1:100], type="o", col=2)

X1 <- sapply(1:q, \(r) y1[r:(n+r)])
X2 <- sapply(1:q, \(r) y2[r:(n+r)])
XX <- cbind(X1,X2)
XXX <- do.call("cbind", lapply(1:q, \(r) cbind(X1[,r],X2[,r])))

res1 <- sapply(1:q, \(k){
  print(k)
  system.time(dac_serial(XX[,c(1:k, q+1:k)], 25))[[3]]
})
res2 <- sapply(1:q, \(k){
  print(k)
  system.time(dac_serial(XXX[,1:(2*k)], 25))[[3]]
})
res3 <- sapply(1:q, \(k){
  print(k)
  # system.time(bruteForceGen(XX[,c(1:k, q+1:k)], serial=T))[[3]]
  system.time(bruteForceGen(XXX[,1:(2*k)], serial=T))[[3]]
})
plot(res1, type="o", ylim = range(c(res1,res2,res3)))
lines(res2, type="o", col=2)
lines(res3, type="o", col=3)



