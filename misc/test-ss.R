library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/dac_serial.cpp")
# sourceCpp("src/dac_serial_debug.cpp")


d <- 30
m <- 2
# n <- 2^18
n <- 300000
X <- replicate(d, 1:n)
dac_serial(X, 100)
system.time(dac_serial(X, 25))[3]



n <- 200
d <- 6
X <- replicate(d, sample(n))
d0 <- sapply(2:d, \(k) bruteForce(X[,1:k]))
d1 <- dac_serial(X, 3)
all(d1 == d0)
# d2 <- dac(X, 3)
# all(d2 == d0)




d <- 2
m <- 2
nn <- m^(10:19)
# n <- 2^18
res <- sapply(nn, \(n){
  print(n)
  mean(replicate(50, {
    X <- replicate(d, sample(n))
    # X <- replicate(d, 1:n)
    system.time(dac_serial(X, 25))[3]
  }))
})
plot(nn, res, type="o")
plot(log(nn,m), log(res,m), type="o")
diff(log(res[c(5,length(res))],m))/diff(log(nn[c(5,length(res))],m))

res2 <- sapply(nn, \(n){
  print(n)
  mean(replicate(50, {
    X <- replicate(d, sample(n))
    # X <- replicate(d, 1:n)
    system.time(pcaPP::cor.fk(X))[3]
  }))
})

plot(nn, res-res[1], type="o")
lines(nn, res2-res2[1], type="o", col=2)
plot(log(nn,m), log(res-res[1]+1,m), type="o")
lines(log(nn,m), log(res2-res2[1]+1,m), type="o", col=2)
diff(log(res[c(5,length(res))],m))/diff(log(nn[c(5,length(res))],m))

nnn <- nn/2^10
nl1 <- nnn*log(nnn)
nl2 <- nnn*log(nnn)^2
plot(nl2, type="o", col=2)
lines(nl1, type="o")

plot(log(nn,m), log(res-res[1]+1,m), type="o")
lines(log(nn,m), log(res2-res2[1]+1,m), type="o", col=2)
lines(log(nn,m), nl1/20000, type="o", col=3)
lines(log(nn,m), nl2/15000, type="o", col=4)
lines(log(nn,m), nl1/2500, type="o", col=3)

