n <- 20
X <- replicate(2, rnorm(n))

th <- pcaPP::cor.fk(X)[1,2]
ths <- sapply(1:nrow(X), \(i) pcaPP::cor.fk(X[-i,])[1,2])

# 1
tts <- n * th - (n-1) * ths
tt <- mean(tts)
sj1 <- sum((tts - tt)^2)/(n-1)
sj2 <- (n-1)*sum((ths - th)^2)

(n-1)*(ths - th) - (tt - tts)

sj1 - sj2

c <- choose(n,2) * (th + 1) / 2
cs <- sapply(ths, \(th) c - choose(n-1,2) * (th + 1) / 2)
sum(cs)
cs - apply(X, 1, \(x) sum(apply(t(X) < x, 2, \(y) identical(y[1],y[2])))-1)

gs <- 2*cs/(n-1) - 1 - th
sh <- 4*sum(gs^2)/n
n*(n-1)*sh/(n-2)^2 - sj2

2*gs/(n-2) - (th - ths)
sh - (n-2)^2*sj2/(n*(n-1))

