library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/dac_serial_simple.cpp")
sourceCpp("src/dac_serial_presort.cpp")

n <- 20
d <- 4
X <- replicate(d, sample(n))
# X <- X[order(X[,1]),]
X
sourceCpp("src/dac_s.cpp")
all(bruteForce(X)==dac_s(X, 3))


dac_serial_P(X, 3)

C0 <- sapply(2:d, \(k) bruteForce(X[,1:k]))
C1 <- dac_serial(X, 3)
C2 <- dac_serial_P(X, 3)

all(C1-C0 == 0)
all(C2-C0 == 0)

sourceCpp("src/dac.cpp")
d1 <- dac(X, 3)
sourceCpp("src/dac.cpp")
d2 <- dac(X, 3)
d0 <- bruteForce(X)
all(d1==d2)

sourceCpp("src/dac_s.cpp")
d4 <- dac_s(X, 3)
all(d1==d2)


sourceCpp("src/dac_serial_simple.cpp")

n <- 10
d <- 4
X <- replicate(d, sample(n))
# X <- replicate(d, 1:n)
# X[,4] <- -X[,4]
system.time(dac_serial(X, 2))
system.time(dac_serial_P(X, 2))
