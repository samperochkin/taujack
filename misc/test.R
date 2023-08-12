library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/test.cpp")

n <- 50
d <- 4
X <- matrix(runif(n*d),n,d)

dac_serial(X,1000)
bruteForce(X)

library(Rcpp)
n <- 30
d <- 4
set.seed(666)
X <- matrix(runif(n*d),n,d)
sourceCpp("src/test4.cpp")
# dac_serial(X,5)
dac_serial(X,5) - dac_serial(X,1000)



sourceCpp("src/bf.cpp")
cbind(bruteForce(X), dac_serial(X,20))



library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/test4.cpp")
sourceCpp("src/dac_serial.cpp")

n <- 30000
d <- 10
set.seed(666)
X <- matrix(runif(n*d),n,d)
# X <- replicate(d, 1:n)
# X[,d] <- -X[,d]
system.time(dac_serial(X,25))
system.time(bruteForce(X))
system.time(pcaPP::cor.fk(X[,1:2]))
system.time(dac_serial(X,100))
