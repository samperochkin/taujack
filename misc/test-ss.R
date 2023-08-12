library(Rcpp)
sourceCpp("src/bf.cpp")


n <- 25
d <- 6
X <- replicate(d, sample(n))

# sourceCpp("src/dst.cpp")
sourceCpp("src/dac_sssserial.cpp")
sourceCpp("src/dac_ssserial.cpp")
sourceCpp("src/dac_sserial.cpp")

d0 <- sapply(2:d, \(k) bruteForce(X[,1:k]))
d1 <- dac_sserial(X, 3)
d2 <- dac_ssserial(X, 3)
d3 <- dac_sssserial(X, 3)
all(d1 == d0)
all(d2 == d0)
all(d3 == d0)

# sourceCpp("src/dac_serial_simple.cpp")
# all(bruteForce(X)==dac_sserial(X, 3)[,d-1])


n <- 75000
d <- 5
X <- replicate(d, sample(n))
system.time(bruteForce(X))
system.time(dac_sserial(X, 25))
system.time(dac_ssserial(X, 25))
system.time(pcaPP::cor.fk(X))


sourceCpp("src/dac_serial_simple.cpp")

n <- 250000
d <- 2
X <- replicate(d, sample(n))
sourceCpp("src/dac_sserial.cpp")
system.time(dac_sserial(X, 25))
sourceCpp("src/dac_serial_simple.cpp")
system.time(dac_sserial(X, 25))


d <- 2
nn <- 1.75^(10:21)
res <- sapply(nn, \(n){
  mean(replicate(50, {
    X <- replicate(d, sample(n))
    sourceCpp("src/dac_sserial.cpp")
    system.time(dac_sserial(X, 25))[3]
  }))
})

plot(log(nn,2), log(res,2), type="o")
