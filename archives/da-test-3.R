
###########################################################################
# Main script (data analysis of daily mean temperatures) ------------------
###########################################################################

# packages ----------------------------------------------------------------
library(tidyverse)
library(lunar)
library(ggplot2)
library(Rcpp)
sourceCpp("src/bf.cpp") # naive/brute force
sourceCpp("src/ms.cpp") # knight's extended
sourceCpp("src/dac_seq.cpp") # divide-and-conquer
source("functions.R") # R wrapper functions and more



# load and arrange data ---------------------------------------------------
data <- readRDS("app/data/data_pp.rds")
data <- data %>% 
  filter(station_name == "Toronto") %>%
  mutate(date = as.Date(date),
         seas = terrestrial.season(date))
         # pseudo_temp = ifelse(seas == "Autumn", pseudo_temp, NA))
stns_id <- unique(data$station_id)
stns_name <- unique(data$station_name)

# reorder to have Toronto in middle
# stns_id <- stns_id[match(c("Welland", "Toronto", "Ottawa"), stns_name)]
# stns_name <- stns_name[match(c("Welland", "Toronto", "Ottawa"), stns_name)]

# check range of data
data %>% group_by(station_name) %>% 
  filter(!is.na(pseudo_temp)) %>%
  summarize(min_date = min(date), max_date = max(date), n = length(date))

# data <- data %>% mutate(pseudo_temp = ifelse(month %in% 3:10, pseudo_temp, NA))

# create dataset of univariate time series
X <- data %>% ungroup %>%
  select(date, station_name, pseudo_temp) %>%
  pivot_wider(names_from = station_name, values_from = pseudo_temp) %>%
  select(any_of(as.character(stns_name))) %>%
  as.matrix()

# dataset for temporal dependence analysis
L <- 4 # number of lags we consider
p <- L+1 # our datasets will thus have p+1 columns
Y <- sapply(1:p, \(r) X[r:(nrow(X)-p+r)])
rm(X); gc()
N <- nrow(Y)

# record which rows to remove when comes the time
na_rows <- which(apply(is.na(Y), 1, any))
n <- N - length(na_rows)




# Data analysis -----------------------------------------------------------

#######################################
#### Pairwise tau at multiple lags ####
#######################################

# pairwise stuff, for comparison
# K <- 365*10
K <- 150
res_pw <- sapply(2:p, \(r) dac_seq(Y[-na_rows, c(1,r)])) |>
  C2Cna(na_rows = na_rows) |> Cna2tsz(K_serial = K, pw = T)
Ths <- res_pw$Th; Ths
res_pw$Sh
matrixcalc::is.positive.definite(res_pw$Sh)

# check zetas
csna <- function(z){i <- which(!is.na(z));zz <- rep(NA, length(z)); zz[i] <- cumsum(z[i]); zz}
Zh <- res_pw$Zh
Zh <- Zh + aperm(Zh, c(2,1,3))
Zhcs <- 4*apply(Zh, c(1,2), \(z) z[1] + csna(z[-1])) |> aperm(perm = c(2,3,1))
i <- 1
par(mfrow = c(2,5), mar=c(2,2,3,1))
# par(mfrow = c(1,1), mar=c(2,2,3,1))
for(j in 1:min(L,10)) plot(Zh[i,j,], main = paste0("pair: (",i,",",j,")"), type="l")
for(j in 1:min(L,10)) plot(Zhcs[i,j,], main = paste0("pair: (",i,",",j,")"), ylim = range(Zhcs[i,1:min(L,10),], na.rm=T), type="l")
mean(Zh[1,10,-(1:100)])
Sh <- Zh2Sh(res_pw$Zh, K=100)
diag(Sh)
matrixcalc::is.positive.definite(Sh)

# plot
df_pw <- data.frame(city = rep(stns_name, each=L),
                    lag = rep(1:L, times=1), p = p,
                    tau = c(Ths),
                    sigma = sqrt(diag(Sh)/n))

ggplot(df_pw, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + #coord_cartesian(xlim=c(0,30)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() + geom_point()

ggplot(df_pw, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + coord_cartesian(xlim=c(1,3)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() + geom_point()


############################################
#### Multivariate tau for each stations ####
############################################

K <- 150
res_multi <- dac_seq(Y[-na_rows,]) |> C2Cna(na_rows = na_rows) |> Cna2tsz(K_serial = K)
Ths <- res_multi$Th; Ths
res_multi$Sh
matrixcalc::is.positive.definite(res_multi$Sh)

# check zetas
csna <- function(z){i <- which(!is.na(z));zz <- rep(NA, length(z)); zz[i] <- cumsum(z[i]); zz}
Zh <- res_multi$Zh
Zh <- Zh + aperm(Zh, c(2,1,3))
Zhcs <- 4*apply(Zh, c(1,2), \(z) z[1] + csna(z[-1])) |> aperm(perm = c(2,3,1))
i <- 1
par(mfrow = c(2,5), mar=c(2,2,3,1))
# par(mfrow = c(1,1), mar=c(2,2,3,1))
for(j in 1:min(L,10)) plot(Zh[i,j,], main = paste0("pair: (",i,",",j,")"), type="l")
for(j in 1:min(L,10)) plot(Zhcs[i,j,], main = paste0("pair: (",i,",",j,")"), ylim = range(Zhcs[i,1:min(L,10),], na.rm=T), type="l")
mean(Zh[1,10,-(1:10)])
Sh <- Zh2Sh(res_multi$Zh, K=20)
diag(Sh)
matrixcalc::is.positive.definite(Sh[1:20,1:20])
#



# Some interesting plots
df_multi <- data.frame(city = rep(stns_name, each=L),
                       lag = rep(1:L, times=1), p = p,
                       tau = c(Ths),
                       sigma = sqrt(diag(Sh)/n))

ggplot(df_multi, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + coord_cartesian(xlim=c(0,30)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() + geom_point()

ggplot(df_multi, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + coord_cartesian(xlim=c(50,65), ylim=c(0,.00001)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() + geom_point()


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
k1 <- which(stns_name == "Toronto")
# k2 <- which(stns_name == "Welland")
k2 <- which(stns_name == "Ottawa")

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
