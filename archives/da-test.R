
###########################################################################
# Main script (data analysis of daily mean temperatures) ------------------
###########################################################################

# packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(Rcpp)
sourceCpp("src/bf.cpp") # naive/brute force
sourceCpp("src/ms.cpp") # knight's extended
sourceCpp("src/dac_seq.cpp") # divide-and-conquer
source("functions.R") # R wrapper functions and more



# load and arrange data ---------------------------------------------------
data <- readRDS("app/data/data_pp_test.rds")
stns_id <- unique(data$station_id)
stns_name <- unique(data$station_name)

# reorder to have Toronto in middle
# stns_id <- stns_id[match(c("Welland", "Toronto", "Ottawa"), stns_name)]
# stns_name <- stns_name[match(c("Welland", "Toronto", "Ottawa"), stns_name)]

# check range of data
data %>% group_by(station_name) %>% 
  filter(!is.na(pseudo_temp)) %>%
  summarize(min_date = min(date), max_date = max(date), n = length(date))

# data <- data %>% mutate(pseudo_temp = ifelse(month %in% c(12,1,2), pseudo_temp, NA))

# create dataset of univariate time series
X <- data %>% ungroup %>%
  select(date, station_name, pseudo_temp) %>%
  pivot_wider(names_from = station_name, values_from = pseudo_temp) %>%
  select(any_of(as.character(stns_name))) %>%
  as.matrix()

# dataset for temporal dependence analysis
L <- 30 # number of lags we consider
p <- L+1 # our datasets will thus have p+1 columns
q <- length(stns_id) # number of stations
Ys <- lapply(1:q, \(k) sapply(1:p, \(r) X[r:(nrow(X)-p+r),k]))
rm(X); gc()
nY <- nrow(Ys[[1]])

# record which rows to remove when comes the time
na_rows0 <- which(apply(is.na(do.call("cbind",Ys)), 1, any))
n0 <- nY - length(na_rows0)
na_rowss <- lapply(1:q, \(k) which(apply(is.na(Ys[[k]][,1:p]), 1, any)))
ns <- nY - sapply(na_rowss, length)

# function to get specific dataset (relies on global env)
pickY <- function(k){
  if(k == 0) return(do.call("cbind", Ys)[-na_rows0,])
  Ys[[k]][-na_rowss[[k]],] 
}



# Data analysis -----------------------------------------------------------

########################################
#### Pairwise taus between stations ####
########################################

# concordance
cols <- (1:q - 1) * p + 1
C <- taujack_ms(pickY(0)[,cols], returnC = T)

# unlist and sub in NAs
C <- apply(combn(1:q,2), 2, \(ij) C[, ij[1], ij[2]]) |> C2Cna(na_rows=na_rows0)

# compute quantities of interest
K <- 400
list2env(Cna2tsz(C, K, pw=T), envir = environment()) # Th, Sh, and Zh
Th; Sh
Zh <- Zh + aperm(Zh, c(2,1,3))
Zhcs <- 2*apply(Zh, c(1,2), \(z) z[1] + 2*cumsum(z[-1])) |> aperm(perm = c(2,3,1))
par(mfrow = c(2,3), mar=c(2,2,3,1))
for(i in 1:q) for(j in i:q) plot(Zhcs[i,j,], main = paste0("pair: (",i,",",j,")"), type="o", cex=.25)
# test of equality will reject for certain... 
# but it's clear there is a problem with zetas.


#########################################################
#### Pairwise tau at multiple lags for each stations ####
#########################################################

# pairwise stuff, for comparison
res_pw <- lapply(1:q, \(k){
  cat("Work on k=", k,".\n")
  i <- (k-1)*p + 1
  C <- sapply(2:p, \(r) dac_seq(pickY(0)[,c(i,i-1+r)])) |>
    C2Cna(na_rows = na_rows0) |> Cna2tsz(K_serial = 400, pw = T)
})
Ths <- sapply(res_pw, "[[", "Th"); Ths
sapply(lapply(res_pw, "[[", "Sh"), matrixcalc::is.positive.definite)

# check zetas
k <- 1
Zh <- res_pw[[k]]$Zh
Zh <- Zh + aperm(Zh, c(2,1,3))
Zhcs <- 2*apply(Zh, c(1,2), \(z) z[1] + 2*cumsum(z[-1])) |> aperm(perm = c(2,3,1))
i <- 1
par(mfrow = c(2,5), mar=c(2,2,3,1))
for(j in 1:min(L,10)) plot(Zh[i,j,], main = paste0("pair: (",i,",",j,")"), type="l")
for(j in 1:min(L,10)) plot(Zhcs[i,j,], main = paste0("pair: (",i,",",j,")"), type="l")
Shs <- lapply(lapply(res_pw, "[[", "Zh"), Zh2Sh, K=15)
sapply(Shs, \(S) matrixcalc::is.positive.definite(S[1:15,1:15]))

# plot
df_pw <- data.frame(city = rep(stns_name, each=L),
                    lag = rep(1:L, times=q), p = p,
                    tau = c(Ths),
                    sigma = sqrt(c(sapply(Shs, \(Sh) diag(Sh)/n0))))

ggplot(df_pw, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + coord_cartesian(xlim=c(0,30)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() + geom_point()


############################################
#### Multivariate tau for each stations ####
############################################

res_multi <- lapply(1:q, \(k){
  cat("Work on k=", k, ".\n")
  dac_seq(pickY(k)) |> C2Cna(na_rows = na_rows0) |> Cna2tsz(K_serial = 400)
})
Ths <- sapply(res_multi, "[[", "Th")

# check zetas
k <- 1
Zh <- res_multi[[k]]$Zh
Zh <- Zh + aperm(Zh, c(2,1,3))
Zhcs <- 2*apply(Zh, c(1,2), \(z) z[1] + 2*cumsum(z[-1])) |> aperm(perm = c(2,3,1))
i <- 1
par(mfrow = c(2,5), mar=c(2,2,3,1))
for(j in 1:min(L,10)){ plot(Zh[i,j,], main = paste0("pair: (",i,",",j,")"), type="l"); abline(h=0, lty=2, col=2)} 
for(j in 1:min(L,10)){ plot(Zh[i,j,], main = paste0("pair: (",i,",",j,")"), type="l", xlim = c(0,50)); abline(h=0, lty=2, col=2)} 
for(j in 1:min(L,10)) plot(Zhcs[i,j,], main = paste0("pair: (",i,",",j,")"), type="l")
Shs <- lapply(lapply(res_multi, "[[", "Zh"), Zh2Sh, K=20)
sapply(Shs, \(S) matrixcalc::is.positive.definite(S[1:15,1:15]))

# Some interesting plots
df_multi <- data.frame(city = rep(stns_name, each=L),
                       lag = rep(1:L, times=q), p = p,
                       tau = c(Ths),
                       sigma = sqrt(c(sapply(Shs, \(Sh) diag(Sh)/ns[k]))))

ggplot(df_multi, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + coord_cartesian(xlim=c(0,30)) +
  geom_line() + geom_point()

ggplot(df_multi, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + coord_cartesian(xlim=c(0,30)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() + geom_point()

ggplot(df_multi, aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_light() + coord_cartesian(xlim=c(20,30), ylim=c(0,.01)) +
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
