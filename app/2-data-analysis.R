
# plot(diff(Ys[[1]][46355 + (-500:500),1]))

# packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(Rcpp)
sourceCpp("src/bf.cpp")
sourceCpp("src/ms.cpp")
sourceCpp("src/dac_serial.cpp")
source("functions.R")
par(mfrow=c(1,1), mar=c(2,2,1,1))

# load and arrange data ---------------------------------------------------
data <- readRDS("app/data/data_pp.rds")
# data <- readRDS("app/data/data_pp_full.rds")
# data <- readRDS("app/data/data_pp4.rds") # ***
# data <- readRDS("app/data/data_pp5.rds") # **
# data <- readRDS("app/data/data_pp9.rds")
# data <- data %>% mutate(seas = lunar::terrestrial.season(date))
stns_id <- unique(data$station_id)[c(2,3,1)]
stns_name <- unique(data$station_name)

X <- data %>% ungroup %>%
  select(date, station_name, pseudo_temp) %>%
  pivot_wider(names_from = station_name, values_from = pseudo_temp) %>%
  select(any_of(as.character(stns_name))) %>%
  as.matrix()
cities <- colnames(X)

# actual dataset
p <- 60 # number of lags we consider
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

# Quick comparison of runtimes --------------------------------------------
runtimes <- function(Y, r){
  t1 <- system.time(dac_serial(Y[,1:r], 25))[3] |> unname()
  t2 <- system.time(bruteForce(Y[,1:r], serial=T))[3] |> unname()
  c(dac = t1, brute_force = t2, ratio = round(t2/t1,2))
}

# in low dimensions the gain is substantial
runtimes(pickY(3), r=2)

# in moderate dimensions the gain is still good, but less impressive
runtimes(pickY(3), r=5)

# for larger dimensions, it seems to be stabilizing
runtimes(pickY(3), r=20)
runtimes(pickY(3), r=30)

# More complete analysis - recomputes all taus everytime. Quite long given we run BF.
# rtimes <- sapply(2:pp, \(k) runtimes(Y[-na_rows2,-(1:pp)], r=k))
# saveRDS(rtimes, "app/objects/rtimes.rds")
rtimes <- readRDS("app/objects/rtimes.rds")
rtimes_df <- data.frame(algorithm = rep(c("divide-and-conquer (DAC)", "brute force (BF)"), times=p),
                        p = rep(2:pp, each=2), runtime = c(rtimes[1:2,]))
ggplot(rtimes_df, aes(x=p-1, y=runtime, shape=algorithm)) + theme_light() +
  geom_line() + geom_point(size=2)

ratio_df <- rtimes_df %>% group_by(p) %>% summarise(ratio = exp(diff(log(runtime))))
print(ratio_df, n = nrow(ratio_df))
diff(rtimes_df[rtimes_df$p == 31,]$runtime)


# Data analysis -----------------------------------------------------------

########################################
#### Pairwise taus between stations ####
########################################

# concordance
cols <- (1:q - 1) * pp + 1
C_temp <- taujack_ms(pickY(0)[,cols], returnC = T)
C0 <- matrix(NA, nY, choose(q,2))
C0[-na_rows0,] <- apply(combn(q,2), 2, \(ij) C_temp[,ij[1],ij[2]])
rm(C_temp)

# tau hat
th0 <- 2*colSums(C0, na.rm=T)/(n0*(n0-1)) - 1
th0

# Hajek (via Sen's Theorem 1, I use K terms) and sigma
g0 <- t((2 * t(C0)/(n0-1) - 1) - th0)
K <- 10
zeta0 <- sapply(0:K, \(k){ 
  gg1 <- g0[1:(nY-k),]; gg2 <- g0[(k+1):nY,]
  sapply(1:nrow(gg1), \(i) tcrossprod(gg1[i,],gg2[i,]), simplify = "array") %>%
    apply(c(1,2), mean, na.rm=T) %>% t()
}, simplify = "array")                
j <- 1
plot(0:K, zeta0[1,j,], type="o", xlim=c(0,K), ylim = range(c(zeta0)))
for(i in 2:q) lines(0:K, zeta0[i,j,], type="o", col=i)
abline(h=0, lty=2)
plot(1:K, cumsum(zeta0[1,1,-1]), type="o", xlim=c(0,K), ylim=range(apply(zeta0, c(2,3), cumsum)))
for(i in 2:q) lines(1:K, cumsum(zeta0[i,i,-1]), type="o", col=i)
abline(h=0, lty=2)
zeta0_sum <- apply(zeta0[,,-1], c(1,2), sum)
Sig0 <- 4*(zeta0[,,1] + zeta0_sum + t(zeta0_sum))
round(Sig0,2)

# quick test of equality will reject for certain... 
# loss <- mahalanobis(sqrt(n0)*(th0[1]-th0[2]), 0, Sig0[1,1]+Sig0[2,2]-2*Sig0[1,2])
# pchisq(loss, 1, lower.tail = F) # barely rejected, perhaps just bad luck.



############################################
#### Multivariate tau for each stations ####
############################################

# pairwise stuff, for comparison
c0m <- colSums(C0, na.rm=T)/(n0*(n0-1))
th0 <- 2*c0m - 1 # as before...
C0s <- lapply(1:q, \(k){
  C_temp <- matrix(NA, nY, p)
  i <- (k-1)*pp + 1
  C_temp[-na_rows0,] <- sapply(2:pp, \(r) dac_serial(pickY(0)[,c(i,i-1+r)]))
  C_temp
})
cm0s <- sapply(1:q, \(k) colSums(C0s[[k]], na.rm=T)/(n0*(n0-1)))
th0s <- 2*cm0s - 1
plot(1:p, th0s[,1], ylim=range(th0s), type="o", cex=.5, pch=19)
for(k in 2:q) lines(1:p, th0s[,k], col=k, type="o", cex=.5, pch=19)
abline(h=0, lty=2)
# g0s <- lapply(1:q, \(k) t(2 * t(C0s[[k]])/(ns[k]-1) - 1 - th0s[,k]))
g0s <- lapply(1:q, \(k) t(2 * t(C0s[[k]])/(n0-1) - 1 - th0s[,k]))
K <- 15
zeta0s <- lapply(1:q, \(i) sapply(0:K, \(k) colMeans(g0s[[i]][1:(nY-k),]*g0s[[i]][(k+1):nY,], na.rm=T)) |> t())
for(k in 1:q){
  plot(0:K, zeta0s[[k]][,1], type="o", xlim=c(0,K))
  for(i in 2:p) lines(0:K, zeta0s[[k]][,i], type="o", col=i)
  abline(h=0, lty=2)
}
csz <- sapply(1:q, \(k) zeta0s[[k]][1,1] + 2*cumsum(zeta0s[[k]][-1,1]))
plot(csz[,1], ylim=range(csz), type="o"); for(k in 2:q) lines(csz[,k], type="o", col=k)
K0 <- 7 # keep only (somewhat arbitrary)
sig0s <- sapply(1:q, \(k) 4*(zeta0s[[k]][1,] + 2*colSums(zeta0s[[k]][1+1:K0,])))
plot(sig0s[,1], ylim=c(0,max(sig0s)), type="o"); for(k in 2:q) lines(sig0s[,k], type="o", col=k)
abline(h=0, lty=2)
sqrt(sig0s[,1])
sqrt(n0)*th0s
plot(pchisq(n0*th0s[,1]^2/sig0s[,1], df=1, lower.tail = F))


# concordance matrices - keep dims of Y to ensure equally spaced timestamps
Cs <- lapply(1:q, \(k){
  C_temp <- matrix(NA, nY, p)
  C_temp[-na_rowss[[k]],] <- dac_serial(pickY(k))
  C_temp
})
cms <- sapply(1:q, \(k) colSums(Cs[[k]], na.rm=T)/(ns[k]*(ns[k]-1)))
tts <- (2^(1:p) * cms - 1)/(2^(1:p) - 1)
gs <- lapply(1:q, \(k) t((2^(1:p) * t(Cs[[k]])/(ns[k]-1) - 1)/(2^(1:p) - 1) - tts[,k]))
plot(1:p, tts[,1], ylim=c(0,max(tts)), type="o", col=1, lwd=1.5)
sapply(2:q, \(k) lines(1:p, tts[,k], col=k, type="o", lty=1, lwd=1.5, pch=k))
sapply(1:q, \(k) lines(1:p, th0s[,k], col=k, type="o", lty=3, lwd=1.5))
abline(h=0, lty=2)

# somehow easier to interpret directly in terms of concordance probability
plot(1:p, cms[,1], ylim=c(0,max(cms)), type="o", col=1, lwd=1.5)
sapply(2:q, \(k) lines(1:p, cms[,k], col=k, type="o", lty=1, lwd=1.5, pch=k))
sapply(1:q, \(k) lines(1:p, cm0s[,k], col=k, type="o", lty=3, lwd=1.5))
abline(h=.5, lty=2)
abline(h=0, lty=2)

# compute the K first terms of the sum in Sen's Theorem 1
K <- 12
zetas <- lapply(1:q, \(i) sapply(0:K, \(k) colMeans(gs[[i]][1:(nY-k),]*gs[[i]][(k+1):nY,], na.rm=T)) |> t())
for(k in 1:q){
  plot(0:K, zetas[[k]][,1], type="o", xlim=c(0,K))
  for(i in 2:p) lines(0:K, zetas[[k]][,i], type="o", col=i)
}

# corresponding jackknife variance
sigs <- sapply(1:q, \(k) 4*(zetas[[k]][1,] + 2*colSums(zetas[[k]][-1,])))
plot(sigs[,1], ylim=c(0,max(sigs)), type="o")
for(k in 2:q) lines(sigs[,k], type="o", col=k)
abline(h=0, lty=2)


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

###############
#### SETUP ####
###############
# Let's focus on autocorrelations with reasonably high variance
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
i <- 1; ii <- (i-1)*(q0-1)+1; j <- 1; jj <- (j-1)*(q0-1)+1
plot(sapply(zetas_test, \(z) z[ii,jj]), ylim=c(0, max(unlist(zetas_test))))
abline(h=0, lty=2)
plot(cumsum(sapply(zetas_test, \(z) z[ii,jj])))
zeta_sum <- Reduce("+", zetas_test[-1])

# Asymptotic var-cov
Sigma_asympt <- 4*(zetas_test[[1]] + zeta_sum + t(zeta_sum))
image(t(Sigma_asympt[(q*q0-q):1,]))

#####################
#### ACTUAL TEST ####
#####################
# cities to compare
k1 <- which(cities == "TORONTO")
# k2 <- which(cities == "WELLAND")
k2 <- which(cities == "OTTAWA CDA")

# corresponding quantities
tts1 <- tts[,k1]
tts2 <- tts[,k2]
ind1 <- (k1-1)*(q0-1) + 1:(q0-1)
ind2 <- (k2-1)*(q0-1) + 1:(q0-1)
S1 <- Sigma_asympt[ind1,ind1]
S2 <- Sigma_asympt[ind2,ind2]
S12 <- Sigma_asympt[ind1,ind2]
n1 <- n_k[[1]][(k1-1)*(q0-1)+1,(k1-1)*(q0-1)+1]
n2 <- n_k[[1]][(k2-1)*(q0-1)+1,(k2-1)*(q0-1)+1]
n12 <- n_k[[1]][(k1-1)*(q0-1)+1,(k2-1)*(q0-1)+1]
n <- n1 + n2 - n12

# var-cov for tt1 - tt2
# the weights ensure that the various sample sizes (and the overlap)
# are appropriately accounted for
Sigma <- n * (S1/n1 + S2/n2 - (S12 + t(S12))*n12/(n1*n2))
Sigma <- (Sigma + t(Sigma))/2 # cancel numerical fuzz
image(t(Sigma[(q0-1):1,]))
plot(diag(Sigma))
image(t((Sigma/sqrt(tcrossprod(diag(Sigma))))[(q0-1):1,]))
eigen(Sigma)$val 
matrixcalc::is.positive.definite(Sigma)

# p-value approximation.
rm_first <- T
if(!rm_first){
  rr0s <- r0s; S <- Sigma
}else{
  rr0s <- r0s[-1]; S <- Sigma[-1,-1]
}
dd <- length(rr0s)
test_stat <- mahalanobis(sqrt(n)*(tts1[rr0s] - tts2[rr0s]), rep(0,dd), S)

# for fun
z_rep <- mvtnorm::rmvnorm(50000, rep(0, dd), S)
test_stat_rep <- mahalanobis(z_rep, rep(0,dd), S)
hist(test_stat_rep, probability=T, breaks=100)
abline(v=test_stat, col=2, lty=2) # not even included (which makes sense given the C.I.)
xx <- seq(0,200,.1)
yy <- dchisq(xx, dd)
lines(xx, yy, col=4, lwd=2)

# Only case not rejected: Toronto-Ottawa with the first coefficient removed.
pval <- c(pchisq(test_stat, df=dd, lower.tail = F), mean(test_stat <= test_stat_rep))
pval
