
# ii <- which.max(abs(g[,1]))
# plot(Ys[[1]][ii + (-500:500),1])
# plot(diff(Ys[[1]][ii + (-500:500),1]))

rr <- 5000
plot((ii + (-rr:rr)) / 365,data[data$station_id == stns_id[1],]$mean_temp[ii + (-rr:rr)], pch=19, cex=.25)
lines((ii + (-rr:rr)) / 365,data[data$station_id == stns_id[1],]$mu[ii + (-rr:rr)], col=2, lwd=2)
abline(v=1:200, lty=2)
plot((ii + (-rr:rr)) / 365,data[data$station_id == stns_id[1],]$ctemp[ii + (-rr:rr)], pch=19, cex=.25)
abline(v=1:200, lty=2)
plot((ii + (-rr:rr)) / 365,data[data$station_id == stns_id[1],]$pseudo_temp[ii + (-rr:rr)], pch=19, cex=.25)
abline(v=1:200, lty=2)
plot((ii + (-rr:(rr-1))) / 365,diff(data[data$station_id == stns_id[1],]$pseudo_temp[ii + (-rr:rr)]), pch=19, cex=.25)
abline(v=1:200, lty=2)

plot(g[1:(nrow(g)-100),1],g[101:nrow(g),1])
plot(g[1:(nrow(g)-1),1],g[2:nrow(g),1])


# packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/ms.cpp")
sourceCpp("src/dac_serial.cpp")
source("functions.R")

# load and arrange data ---------------------------------------------------
data <- readRDS("app/data/data_pp.rds")
# data <- data %>% mutate(seas = lunar::terrestrial.season(date))
stns_id <- unique(data$station_id)
stns_name <- unique(data$station_name)

X <- data %>% ungroup %>%
  select(date, station_name, pseudo_temp) %>%
  pivot_wider(names_from = station_name, values_from = pseudo_temp) %>%
  select(any_of(as.character(stns_name))) %>%
  as.matrix()
cities <- colnames(X)

# actual dataset
p <- 30 # number of lags we consider
pp <- p+1 # our datasets will thus have p+1 columns
q <- length(stns_id) # number of stations
Ys <- lapply(1:q, \(k) sapply(1:pp, \(r) X[r:(nrow(X)-pp+r),k]))
rm(X); gc()
nY <- nrow(Ys[[1]])

# record which rows to remove when comes the time
na_rows0 <- which(apply(is.na(do.call("cbind",Ys)), 1, any))
n0 <- nY - length(na_rows0)
na_rowss <- lapply(1:q, \(k) which(apply(is.na(Ys[[k]][,1:pp]), 1, any)))
ns <- nY - sapply(na_rowss, length)

# function to get specific dataset
pickY <- function(k, YYs = Ys, nna_rowss = na_rowss, nna_rows0 = na_rows0){
  if(k == 0) return(do.call("cbind", Ys)[-na_rows0,])
  Ys[[k]][-na_rowss[[k]],] 
}


# Data analysis -----------------------------------------------------------

############################################
#### Multivariate tau for each stations ####
############################################

# concordance matrices - keep dims of Y to ensure equally spaced timestamps
p <- 30
pp <- 31
C <- matrix(NA, nY, p)
C[-na_rowss[[1]],] <- dac_serial(pickY(1))
cm <- colSums(C, na.rm=T)/(ns[1]*(ns[1]-1))
tt <- (2^(1:p) * cm - 1)/(2^(1:p) - 1)
g <- t((2^(1:p) * t(C)/(ns[1]-1) - 1)/(2^(1:p) - 1) - tt)
plot(1:p, tt, ylim=c(0,max(tt)), type="o", col=1, lwd=1.5)


C_jack <- t((2^(1:p) * (colSums(C, na.rm=T) - t(C))/((ns[1]-1)*(ns[1]-2)) - 1)/(2^(1:p) - 1))
hist(C_jack[,1], breaks=100)
dim(C_jack)

# compute the K first terms of the sum in Sen's Theorem 1
K <- 300
zetas <- sapply(0:K, \(k) colMeans(g[1:(nY-k),]*g[(k+1):nY,], na.rm=T)) |> t()
plot(0:K, zetas[,1], type="o", xlim=c(0,K))
for(i in 2:p) lines(0:K, zetas[,i], type="o", col=i)


# K <- 1000
# r <- 20
# # zz <- colMeans(g[1:(nY-K),]*g[(K+1):nY,], na.rm=T)
# zz <- sapply(0:K, \(k) mean(g[1:(nY-k),r]*g[(k+1):nY,r], na.rm=T))
# plot(0:K/365.25, zz)
# abline(v=1:10, lty=2)
# plot(0:K/365.25, cumsum(zz))
# abline(v=1:10, lty=2)
# 
# K0 <- which(zz < 0)[1]
# mean(zz[K0:K])
# zzz <- zz
# zzz[K0:K] <- mean(zz[K0:K])
# zzz[K0:K] <- 1e-7
# plot(cumsum(zzz))
# 
# mean(zz[750:K])

# corresponding jackknife variance
sigs <- 4*(zetas[1,] + 2*colSums(zetas[-1,]))
plot(sigs, ylim=c(0,max(sigs)), type="o")


# Some interesting plots
df <- data.frame(city = rep(cities, each=p),
                 lag = rep(1:p, times=q), p = p,
                 tau = c(tts),
                 sigma = sqrt(c(sapply(1:q, \(k) sigs[,k]/ns[k]))))

ggplot(df, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + coord_cartesian(xlim=c(0,30)) +
  geom_line() + geom_point()

ggplot(df, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + coord_cartesian(xlim=c(0,30)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() + geom_point()

ggplot(df, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + coord_cartesian(xlim=c(1,4), ylim=c(.3,.45)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() + geom_point()

ggplot(df, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() +
  coord_cartesian(xlim=c(1,6), ylim=c(.05,.42)) +
  # coord_cartesian(xlim=c(0,5), ylim=c(.1,.5)) +
  # coord_cartesian(xlim=c(15,30), ylim=c(0,.018)) +
  # coord_cartesian(xlim=c(45,60), ylim=c(0,.0001)) +
  # coord_cartesian(xlim=c(30,60), ylim=c(-.000002,.00001)) +
  # coord_cartesian(xlim=c(4,15), ylim=c(0,.2)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() +
  geom_point()




# Inference (test of equality of the two sets autocorrelations) -----------
# Let's focus on autocorrelations with reasonably high variance
# Those with extremely low variance are likely to pollute the inference
r0s <- 1:15
q0 <- length(r0s)+1
g_test <- do.call("cbind", lapply(1:q, \(k) gs[[k]][, r0s]))
K <- 30 # number of terms in sum (Sen)

# Record various sample sizes to correctly compute the Sigma matrix 
n_k <- lapply(0:K, \(k) t(!is.na(g_test)[1:(nY-k),]) %*% !is.na(g_test)[(k+1):nY,])
sapply(1:q, \(k) ns[k] == n_k[[1]][(k-1)*q0 + 1, (k-1)*q0 + 1])  # check
g_test[is.na(g_test)] <- 0 # can set that to zero since we know the samples sizes

# Computation of the summands for the variance-covariance matrix
zetas_test <- lapply(0:K, \(k) (t(g_test[1:(nY-k),]) %*% g_test[(k+1):nY,])/n_k[[k+1]])
par(mfrow=c(1,1), mar=c(2,2,1,1))
j <- 2
plot(sapply(zetas_test, \(z) z[1,j]), ylim=c(0, max(unlist(zetas_test))))
abline(h=0, lty=2)
plot(cumsum(sapply(zetas_test, \(z) z[1,j])))
zeta_sum <- Reduce("+", zetas_test[-1])
Sigma_asympt <- 4*(zetas_test[[1]] + zeta_sum + t(zeta_sum))
image(t(Sigma_asympt)[(q*q0-q):1,])

# Compute Sigma
k1 <- which(cities == "OTTAWA CDA")
tts1 <- tts[,k1]
ind1 <- (k1-1)*(q0-1) + 1:(q0-1)
S1 <- Sigma_asympt[ind2,ind2]
# n1 <- ns[k1]
n1 <- n_k[[1]][(k1-1)*(q0-1)+1,(k1-1)*(q0-1)+1]

k2 <- which(cities == "TORONTO")
tts2 <- tts[,k2]
ind2 <- (k2-1)*(q0-1) + 1:(q0-1)
S2 <- Sigma_asympt[ind2,ind2]
# n2 <- ns[k2]
n2 <- n_k[[1]][(k2-1)*(q0-1)+1,(k2-1)*(q0-1)+1]

S12 <- Sigma_asympt[ind1,ind2]
n <- min(n1,n2)

# Note the weights to ensure that the various sample sizes (and the overlap) are accounted for
Sigma <- (n1 + n2) * (S1/n1 + S2/n2 - (sqrt(n)/(sqrt(n1)*sqrt(n2)))^2 * (S12 + t(S12)))
# Sigma <- S1 + S2 - S12 - t(S12)
Sigma <- (Sigma + t(Sigma))/2 # cancel numerical fuzz
image(t(Sigma)[(q0-1):1,])
plot(diag(Sigma))
image(t(Sigma/sqrt(tcrossprod(diag(Sigma))))[(q0-1):1,])
eigen(Sigma)$val 
matrixcalc::is.positive.definite(Sigma)
Sigma_inv <- solve(Sigma)

# p-value approximation -- Note the factor (n1 + n2), consistent with that in Sigma.
test_stat <- mahalanobis(sqrt(n1+n2)*(tts1[r0s] - tts2[r0s]), rep(0,q0-1), Sigma_inv, TRUE)

# mc replicates just for fun
z_rep <- mvtnorm::rmvnorm(50000, rep(0, q0-1), Sigma)
test_stat_rep <- mahalanobis(z_rep, rep(0,q0-1), Sigma_inv, TRUE)

hist(test_stat_rep, probability=T, breaks=100)
abline(v=test_stat, col=2, lty=2) # not even included (which makes sense given the C.I.)
xx <- seq(0,200,.1)
yy <- dchisq(xx, q0-1)
lines(xx, yy, col=4, lwd=2)

pval <- c(pchisq(test_stat, df=q0-1, lower.tail = F),
          mean(test_stat <= test_stat_rep))
pval









n <- n_k[[1]][(k1-1)*(q0-1)+1,(k2-1)*(q0-1)+1]

# Note the weights to ensure that the various sample sizes (and the overlap) are accounted for
Sigma <- (n1 + n2 - n) * (S1/n1 + S2/n2 - (S12 + t(S12))/n)
# Sigma <- S1 + S2 - S12 - t(S12)
Sigma <- (Sigma + t(Sigma))/2 # cancel numerical fuzz
image(t(Sigma)[(q0-1):1,])
plot(diag(Sigma))
image(t(Sigma/sqrt(tcrossprod(diag(Sigma))))[(q0-1):1,])
eigen(Sigma)$val 
matrixcalc::is.positive.definite(Sigma)
Sigma_inv <- solve(Sigma)

# p-value approximation -- Note the factor (n1 + n2), consistent with that in Sigma.
test_stat <- mahalanobis(sqrt(n1+n2-n)*(tts1[r0s] - tts2[r0s]), rep(0,q0-1), Sigma_inv, TRUE)

# mc replicates just for fun
z_rep <- mvtnorm::rmvnorm(50000, rep(0, q0-1), Sigma)
test_stat_rep <- mahalanobis(z_rep, rep(0,q0-1), Sigma_inv, TRUE)

hist(test_stat_rep, probability=T, breaks=100)
abline(v=test_stat, col=2, lty=2) # not even included (which makes sense given the C.I.)
xx <- seq(0,200,.1)
yy <- dchisq(xx, q0-1)
lines(xx, yy, col=4, lwd=2)

pval <- c(pchisq(test_stat, df=q0-1, lower.tail = F),
          mean(test_stat <= test_stat_rep))
pval

