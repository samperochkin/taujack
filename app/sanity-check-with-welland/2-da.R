
# packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/dac_serial.cpp")


# load and arrange data ---------------------------------------------------
data <- readRDS("app/data/data_pp_alt.rds")
# data <- data %>% mutate(seas = lunar::terrestrial.season(date))
stns_id <- unique(data$station_id)

X <- data %>% ungroup %>%
  select(date, station_id, pseudo_temp) %>%
  pivot_wider(names_from = station_id, values_from = pseudo_temp) %>%
  select(any_of(as.character(stns_id))) %>%
  as.matrix()

# X[is.na(X[,1]),2] <- NA
# X[is.na(X[,2]),1] <- NA

# Quick comparison of runtimes --------------------------------------------

# window we will look at
d <- ncol(X) # d instead of p, confusing sorry
q <- 30

# dataset - note the order (station=1 lag=1, station=2 lag=1, station=1 lag=2, ...)
Y <- lapply(1:q, \(r){
  sapply(1:d, \(k) X[r:(nrow(X)-q+r), k])
}) %>% do.call(what="cbind")
na_rows <- which(apply(is.na(Y), 1, any))
nY <- nrow(Y)
n <- nY - length(na_rows)

# function to compute runtimes
runtimes <- function(Y, r){
  t1 <- system.time(dac_serial(Y[,1:(2*r)]))[3] |> unname()
  t2 <- system.time(bruteForce(Y[,1:(2*r)], serial=T))[3] |> unname()
  c(dac = t1, brute_force = t2, ratio = round(t2/t1,2))
}

# in low dimensions the gain is substantial
runtimes(Y[-na_rows,], r=1)

# in moderate dimensions the gain is still good, but less impressive
runtimes(Y[-na_rows,], r=5)

# for larger dimensions, it seems to be stabilizing
runtimes(Y[-na_rows,], r=20)
runtimes(Y[-na_rows,], r=30)




# Data analysis -----------------------------------------------------------
# recall that (or note the change)
d <- ncol(X)
q <- 5

# dataset - note the order (station=1 lag=1, station=2 lag=1, station=1 lag=2, ...)
# This is to allow computation of tau_tilde (in the paper) when p>1
Y <- lapply(1:q, \(r){
  sapply(1:d, \(k) X[r:(nrow(X)-q+r), k])
}) %>% do.call(what="cbind")
na_rows <- which(apply(is.na(Y), 1, any))
nY <- nrow(Y)
n <- nY - length(na_rows)
na_rows1 <- which(apply(is.na(Y[,2*(1:q)-1]), 1, any))
n1 <- nY - length(na_rows1)
na_rows2 <- which(apply(is.na(Y[,2*(1:q)]), 1, any))
n2 <- nY - length(na_rows2)

# joint concordance - double the columns and then aggregate over spatial dim.
C <- matrix(NA, nY, q)
C[-na_rows,] <- dac_serial(Y[-na_rows,])[,seq(1,d*q,2)]
kks <- seq(2,2*q,2)
ths <- (2^(kks-1) * colSums(C, na.rm=T)/(n*(n-1)) - 1)/(2^(kks-1) - 1)
g <- t((2^(kks-1) * t(C)/(n-1) - 1)/(2^(kks-1) - 1) - ths)
K <- 10 # Set this high, and then adjust according to later plots
zetas <- sapply(0:K, \(k) colMeans(g[1:(nY-k),]*g[(k+1):nY,], na.rm=T)) |> t()
sigs <- zetas[1,] + 2*colSums(zetas[-1,])

# univariate concordance
C1 <- C2 <- matrix(NA, nY, q-1)
C1[-na_rows1,] <- dac_serial(Y[-na_rows1,2*(1:q)-1])
C2[-na_rows2,] <- dac_serial(Y[-na_rows2,2*(1:q)])
ths1 <- (2^(2:q-1) * colSums(C1, na.rm=T)/(n1*(n1-1)) - 1)/(2^(2:q-1) - 1)
ths2 <- (2^(2:q-1) * colSums(C2, na.rm=T)/(n2*(n2-1)) - 1)/(2^(2:q-1) - 1)
g1 <- t((2^(2:q-1) * t(C1)/(n1-1) - 1)/(2^(2:q-1) - 1) - ths1)
g2 <- t((2^(2:q-1) * t(C2)/(n2-1) - 1)/(2^(2:q-1) - 1) - ths2)
zetas1 <- sapply(0:K, \(k) colMeans(g1[1:(nY-k),]*g1[(k+1):nY,], na.rm=T)) |> t()
zetas2 <- sapply(0:K, \(k) colMeans(g2[1:(nY-k),]*g2[(k+1):nY,], na.rm=T)) |> t()
sigs1 <- zetas1[1,] + 2*colSums(zetas1[-1,])
sigs2 <- zetas2[1,] + 2*colSums(zetas2[-1,])

# check taus
plot(0:(q-1), ths, ylim=c(0,max(c(ths,ths1,ths2))), type="o")
lines(1:(q-1), ths1, col=2, type="o")
lines(1:(q-1), ths2, col=3, type="o")
abline(h=0, lty=2)

# check zetas and sigs
plot(0:K, zetas[,1], type="o", xlim=c(0,K))
for(i in 2:(q-1)) lines(0:K, zetas[,i], type="o", col=i)
abline(h=0, lty=2)
plot(sigs, ylim=c(0,max(sigs)), type="o")
abline(h=0, lty=2)

plot(0:K, zetas1[,1], type="o", xlim=c(0,K))
for(i in 2:(q-1)) lines(0:K, zetas1[,i], type="o", col=i)
abline(h=0, lty=2)
plot(sigs1, ylim=c(0,max(sigs1)), type="o")
abline(h=0, lty=2)

plot(0:K, zetas2[,1], type="o", xlim=c(0,K))
for(i in 2:(q-1)) lines(0:K, zetas2[,i], type="o", col=i)
abline(h=0, lty=2)
plot(sigs2, ylim=c(0,max(sigs2)), type="o")
abline(h=0, lty=2)


# plots
cities <- c("ctiy1", "ctiy2")
df <- data.frame(city = c("joint", rep(c("joint", cities), each=q-1)),
                 lag = c(0:(q-1),1:(q-1),1:(q-1)),
                 tau = c(ths, ths1, ths2),
                 sigma = sqrt(c(sigs/n,sigs1/n1,sigs2/n2)),
                 p = d) # yeah, a bit confusing sorry.

ggplot(df, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() +
  # coord_cartesian(xlim=c(8,20), ylim=c(0,.1)) +
  # geom_ribbon(aes(ymin = tau - 3*sigma, ymax = tau + 3*sigma), alpha=.3, size=0) +
  geom_line() +
  geom_point()

ggplot(df, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() +
  # coord_cartesian(xlim=c(8,20), ylim=c(0,.1)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() +
  geom_point()

ggplot(df, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() +
  # coord_cartesian(xlim=c(15,30), ylim=c(0,.018)) +
  # coord_cartesian(xlim=c(45,60), ylim=c(0,.0001)) +
  coord_cartesian(xlim=c(30,60), ylim=c(-.000002,.00001)) +
  geom_ribbon(aes(ymin = tau - 3*sigma, ymax = tau + 3*sigma), alpha=.3, size=0) +
  geom_line() +
  geom_point()



# Inference (test of equality of the two sets autocorrelations) -----------
g_test <- cbind(g1,g2)

# Record various sample sizes to correctly compute the Sigma matrix 
n_k <- lapply(0:K, \(k){
  t(!is.na(g_test)[1:(nY-k),]) %*% !is.na(g_test)[(k+1):nY,]
})
n == n_k[[1]][1,q] # check
n1 == n_k[[1]][1,1] # check
n2 == n_k[[1]][q,q] # check
g_test[is.na(g_test)] <- 0 # can set that to zero since we know the counts

# Computation of the summands for the variance-covariance matrix
zetas_test <- lapply(0:K, \(k){
  # M <- (t(g_test[1:(nY-k),]) %*% g_test[(k+1):nY,]) # we don't divide yet
  # (M + t(M))/max(n1,n2) # to conserve approppriate ratios
  M <- (t(g_test[1:(nY-k),]) %*% g_test[(k+1):nY,])#/n_k[[k+1]] # we don't divide yet
  (M + t(M))/(n_k[[k+1]]+t(n_k[[k+1]])) # corresponds to weighted average
  # M <- (t(g_test[1:(nY-k),]) %*% g_test[(k+1):nY,])/n_k[[k+1]]
  # (M + t(M))
})
# Note the absence of the factor two in front of Reduce below.
# It is already accounted for in M + t(M) just above.
Sigma_asympt <- 4*(zetas_test[[1]] + 2*Reduce("+", zetas_test[-1]))
image(t(Sigma_asympt)[(2*q-2):1,])

# Sigma_asympt2 <- min(n1,n2)*Sigma_asympt/n_k[[1]]
Sigma_asympt2 <- n_k[[1]]*Sigma_asympt/(n1 + n2)
image(t(Sigma_asympt2)[(2*q-2):1,])

# Compute pooled estimate of Sigma
ind1 <- 1:(q-1); ind2 <- q-1 + 1:(q-1)
S1 <- Sigma_asympt2[ind1,ind1]
S2 <- Sigma_asympt2[ind2,ind2]
S12 <- Sigma_asympt2[ind1,ind2]
# Sigma <- ((n1+n2)/2)*(S1/n1 + S2/n1 - (S12 + t(S12))/n)
Sigma <- S1 + S2 - S12 - t(S12)
Sigma <- (Sigma + t(Sigma))/2 # cancel numerical fuzz
image(t(Sigma)[(q-1):1,])
eigen(Sigma)$val 
matrixcalc::is.positive.definite(Sigma)
Sigma_inv <- solve(Sigma)

# p-value approximation -- Monte Carlo, but this is probably a chi squared or very close.
test_stat <- mahalanobis(sqrt(n1+n2)*(ths1 - ths2), rep(0,q-1), Sigma_inv, TRUE)
z_rep <- mvtnorm::rmvnorm(50000, rep(0, q-1), Sigma)
test_stat_rep <- mahalanobis(z_rep, rep(0,q-1), Sigma_inv, TRUE)

hist(test_stat_rep, probability=T, breaks=100)
abline(v=test_stat, col=2, lty=2) # not even close (which makes sense given the C.I.)
xx <- seq(0,200,.1)
yy <- dchisq(xx, q-1) # It's probably what this is...
lines(xx, yy, col=4, lwd=2)

pval <- pchisq(test_stat, df=q-1, lower.tail = F)
pval # pretty much zero.
