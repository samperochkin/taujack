
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
data <- readRDS("app/data/data_pp.rds")
stns_id <- unique(data$station_id)
stns_name <- unique(data$station_name)

# reorder to have Toronto in middle
stns_id <- stns_id[match(c("Welland", "Toronto", "Ottawa"), stns_name)]
stns_name <- stns_name[match(c("Welland", "Toronto", "Ottawa"), stns_name)]

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
L <- 29 # number of lags we consider
p <- L+1 # our datasets will thus have p+1 columns
q <- length(stns_id) # number of stations
Ys <- lapply(1:q, \(k) sapply(1:p, \(r) X[r:(nrow(X)-p+r),k]))
rm(X); gc()
N <- nrow(Ys[[1]])

# record which rows to remove when comes the time
na_rows0 <- which(apply(is.na(do.call("cbind",Ys)), 1, any))
n0 <- N - length(na_rows0)
na_rowss <- lapply(1:q, \(k) which(apply(is.na(Ys[[k]][,1:p]), 1, any)))
ns <- N - sapply(na_rowss, length)

# function to get specific dataset (relies on global env)
pickY <- function(k){
  if(k == 0) return(do.call("cbind", Ys)[-na_rows0,])
  Ys[[k]][-na_rowss[[k]],] 
}



# Quick comparison of runtimes (w/ Toronto) -------------------------------
runtimes <- function(Y, r){
  cat("computing tau for dataset with n=", nrow(Y), "obs. and p=", ncol(Y), "dimensions, using", r, "of them.\n")
  t1 <- system.time(dac_seq(Y[,1:r], 25))[3] |> unname()
  t2 <- system.time(bruteForce(Y[,1:r], seq=T))[3] |> unname()
  c(dac = t1, brute_force = t2, ratio = round(t2/t1,2))
}

# in low dimensions the gain is substantial
runtimes(pickY(2), r=2)

# in moderate dimensions the gain is still very good, but less impressive
runtimes(pickY(2), r=5)

# for larger dimensions, it seems to be stabilizing
runtimes(pickY(2), r=20)
runtimes(pickY(2), r=30)

# More complete analysis - recomputes all taus everytime. Quite long given we run BF.
# rtimes <- sapply(2:p, \(k) runtimes(pickY(2), r=k))
# saveRDS(rtimes, "app/rtimes.rds")
rtimes <- readRDS("app/rtimes.rds")
rtimes_df <- data.frame(algorithm = rep(c("DAC", "BF"), times=p-1),
                        p = rep(2:p, each=2), runtime = c(rtimes[1:2,]))
ggplot(rtimes_df, aes(x=p-1, y=runtime, linetype=algorithm)) +
  theme_bw() + xlab("maximum lag (p-1)") + ylab("runtime (in sec.)") +
  scale_linetype_manual(values = c(1,2)) +
  geom_line(size=.25)# + geom_point(size=.25)
ggsave("figures/runtimes.pdf", device = "pdf", width = 6.5, height = 2, units = "in")

ratio_df <- rtimes_df %>% group_by(p) %>% summarise(ratio = exp(diff(log(runtime))))
print(ratio_df, n = nrow(ratio_df))
diff(rtimes_df[rtimes_df$p == 30,]$runtime)



# Data analysis -----------------------------------------------------------

# C2Cna:  Function that takes the concordance matrix C and embeds it
#         in a matrix that contains NAs representing missing obs.
#         This is ease the computation of the zetas, which requires
#         knowing the timestamps associated with the counts

# Cna2tsz:  Computes the quantities of interest from C (w/ NAs).
#           It returns Th (the taus), Sh (jack var-cov) and Zh (Z, zetas)

# Zh2Sh: Constructs Sh out of Zh, using K (given as input) terms



#########################################################
#### Pairwise tau at multiple lags for each stations ####
#########################################################

# omitted from the analysis in the paper.

K <- 100 # number of terms of the zeta sum we want to compute
res <- lapply(1:q, \(k){ # Cna2tsz takes a bit of time, but nothing crazy
  cat("Work on k=", k,".\n")
  i <- (k-1)*p + 1 # desired time series (r, below, is the lagged one)
  sapply(2:p, \(r) dac_seq(pickY(0)[,c(i,i-1+r)])) |>
    C2Cna(na_rows = na_rows0) |> Cna2tsz(K_serial = K, pw = T)
})
Ths <- sapply(res, "[[", "Th")
Shs <- lapply(res, "[[", "Sh")
Zhs <- lapply(res, "[[", "Zh")
rm(res)

# quick look
Ths
sapply(Shs, matrixcalc::is.positive.definite)

####---- select number of terms we want to use in the zeta sum 
k <- 1 # station (change as wished, from 1 to 3)
i <- 1 # we focus on tau_ij (in the plots only)
par(mfrow = c(2,5), mar=c(2,2,3,1))
Zh <- Zhs[[k]]
for(j in 1:min(L,10)){
  plot(Zh[i,j,], main = paste0("pair: (",i,",",j,")"), type="l", ylim = range(Zh[i,1:min(L,10),]))
  abline(h=0, lty=2, col=2)
} 

# computes the series of partial sums of Z + t(Z), to which Z0 is added
Shcs <- 4*apply(Zh, c(1,2), \(z) z[1] + cumsum(z[-1])) |> aperm(perm = c(2,3,1))
for(j in 1:min(L,10)) plot(Shcs[i,j,], main = paste0("pair: (",i,",",j,")"),
                           ylim = range(Shcs[i,1:min(L,10),], na.rm=T), type="l")
# these plots look like what I have seen on synthetic data
par(mfrow = c(1,1), mar=c(2,2,3,1))
####-----------------------------

# update Shs so that it uses K=20 terms only
Shs <- lapply(Zhs, Zh2Sh, K=20)
sapply(Shs, matrixcalc::is.positive.definite)

# plot (with choice of zoom)
cc <- NULL # no zoom
# cc <- coord_cartesian(xlim = c(1,1.5), ylim = c(.32,.49)) # lag-1
# cc <- coord_cartesian(xlim = c(2,5), ylim = c(.07,.3)) # lag-2 to 5
# cc <- coord_cartesian(xlim = c(L-5,L), ylim = c(0,.04)) # tail

data.frame(city = rep(stns_name, each=L), lag = rep(1:L, times=q), p = p,
           tau = c(Ths), sigma = sqrt(c(sapply(Shs, \(Sh) diag(Sh)/n0)))) %>%
  ggplot(aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_bw() + cc +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, size=0) +
  geom_line() + geom_point()



############################################
#### Multivariate tau for each stations ####
############################################

K <- 100 # number of terms of the zeta sum we want to compute
res <- lapply(1:q, \(k){ # Cna2tsz takes a bit of time, but nothing crazy
  cat("Work on k=", k,".\n")
  dac_seq(pickY(k)) |> C2Cna(na_rows = na_rowss[[k]]) |> Cna2tsz(K_serial = K)
})
Ths <- sapply(res, "[[", "Th")
Shs <- lapply(res, "[[", "Sh")
Zhs <- lapply(res, "[[", "Zh")

# quick look
Ths
sapply(Shs, matrixcalc::is.positive.definite)
# when we focus on the first ten components to avoid numerical instability
sapply(Shs, \(S) matrixcalc::is.positive.definite(S[1:10, 1:10]))

####---- select number of terms we want to use in the zeta sum 
k <- 1 # station (change as wished)
i <- p # we focus on tau_ij (in the plots only)
par(mfrow = c(2,5), mar=c(2,2,3,1))
Zh <- Zhs[[k]]
for(j in 1:min(L,10)){
  plot(Zh[i,j,], main = paste0("pair: (",i,",",j,")"), type="l")
  abline(h=0, lty=2, col=2)
} 

# computes the series of partial sums of Z, to which Z0 is added
Shcs <- 4*apply(Zh, c(1,2), \(z) z[1] + cumsum(z[-1])) |> aperm(perm = c(2,3,1))
for(j in 1:min(L,10)){
  plot(Shcs[i,j,], main = paste0("pair: (",i,",",j,")"), ylim = range(Shcs[i,1:min(L,10),], na.rm=T), type="l")
  abline(h=0, lty=2, col=2)
}

# looks like similar plots from synthetic data
par(mfrow = c(1,1), mar=c(2,2,3,1))
####-----------------------------

# update Shs so that it uses K=20 terms only
# here we focus on the first ten components to avoid numerical instability
Shs <- lapply(Zhs, Zh2Sh, K=20)
sapply(Shs, \(S) matrixcalc::is.positive.definite(S[1:10, 1:10]))

# plot (below is choice of zoom)
cc <- NULL # none
# cc <- coord_cartesian(xlim = c(1,2.5), ylim = c(.35,.49)) # lag-1-2
# cc <- coord_cartesian(xlim = c(2,5), ylim = c(.18,.45)) # lag-2 to 5
# cc <- coord_cartesian(xlim = c(L-5,L), ylim = c(0,.004)) # tail

data.frame(city = rep(stns_name, each=L), lag = rep(1:L, times=q), p = p,
           tau = c(Ths),
           sigma = sqrt(c(sapply(1:q, \(k) diag(Shs[[k]])/ns[k])))) %>%
  ggplot(aes(x = lag, y = tau, col=city, group=city, fill=city)) +
  theme_bw() + cc +
  xlab("lag (k-1)") + ylab(bquote(hat(tau)["k"])) +
  theme(legend.position = c(.75,.6),
        text = element_text(size=9),
        axis.text=element_text(size=9),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.background = element_rect(colour="black"),
        legend.title = element_text(size = 9),
        legend.text=element_text(size=9)) +
  geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, colour=NA) +
  geom_line(size=.25) + geom_point(size=.25)
# for figure in paper, see tau-figure.R


# Inference (test of equality of the two sets autocorrelations) -----------
# In what follows, we mimic what Cna2tsz does, but we account for
# the overlap between some of the data from distinct time series

# Let's focus on the first 10 lags, for convenience
L <- 15
p <- L+1
K <- 20

# collect gh from our res object (computed earlier)
gh <- do.call("cbind", lapply(1:q, \(k) res[[k]]$gh[,1:L]))
Th <- colMeans(gh, na.rm=T)

# nn0[i,j] is the total number of samples for Sh[i,j]
nn0 <- t(!is.na(gh)) %*% !is.na(gh)

Zh <- sapply(0:K, \(j){
  if(j %% 5 == 0) cat(round(100*j/(K + 1),1), "% -- ")
  
  G1 <- gh[1:(N-j),,drop=F]
  G2 <- gh[(j+1):N,,drop=F]
  
  # This is an important part
  # It records the number of joint samples for that specific lag value
  nn <- t(!is.na(G1)) %*% !is.na(G2)
  if(any(nn == 0)) stop("There's something to be taken care of...")
  
  # weight to take into account that cov divides by n-1, rather than n (which we want)
  w1 <- (nn-1)/nn
  
  # weight for the terms in the zeta sum, so that the contribution better
  # match in the empirical case
  w2 <- nn/nn0
  
  # this is to get the asymptotic Sigma matrix:
  z <- w1 * w2 * cov(G1, G2, use = "pairwise") # we use all we can here
  
  # "twice the covariance" (not for Z_0)
  if(k>0) z <- z + t(z)
  
  return(z)
}, simplify = "array")

Zsum <- apply(Zh[,,-1], c(1,2), sum, na.rm=T)
Sh <- 4*(Zh[,,1] + Zsum)
image(t(Sh[(q*p-q):1,]))
matrixcalc::is.positive.definite(Sh)

#####################
#### ACTUAL TEST ####
#####################
stns_name

# collect index
ind1 <- (1-1)*(p-1) + 1:(p-1)
ind2 <- (2-1)*(p-1) + 1:(p-1)
ind3 <- (3-1)*(p-1) + 1:(p-1)

# keep only...
keep <- 1 # first component (p=2, first test in paper)
# keep <- 1:2 # first two component (p=3, second test in paper)
# keep <- 1:11 # first 11 components (p=12, third test in paper)
ind1 <- ind1[keep]; ind2 <- ind2[keep]; ind3 <- ind3[keep]
l <- length(keep)

# submatrices of Sigma with corresponding samples sizes
S1 <- Sh[ind1,ind1]; n1 <- nn0[ind1[1],ind1[1]]
S2 <- Sh[ind2,ind2]; n2 <- nn0[ind2[1],ind2[1]]
S3 <- Sh[ind3,ind3]; n3 <- nn0[ind3[1],ind3[1]]
S12 <- Sh[ind1,ind2]; n12 <- nn0[ind1[1],ind2[1]]
S23 <- Sh[ind2,ind3]; n23 <- nn0[ind2[1],ind3[1]]

# compute new Sigma matrix adjusting for sample sizes
# based on the expansion of sqrt(N)*(thi - thj)
Sig12 <- N*(S1/n1 + S2/n2 - (S12 + t(S12))*n12/(n1*n2))
Sig23 <- N*(S2/n2 + S3/n3 - (S23 + t(S23))*n23/(n2*n3))
matrixcalc::is.positive.definite(Sig12)
matrixcalc::is.positive.definite(Sig23)
image(t(Sig12[l:1,]))
image(t(Sig23[l:1,]))

# test stat and test (w/ Monte Carlo too, for fun)
test_stat12 <- mahalanobis(sqrt(N)*(Th[ind1] - Th[ind2]), rep(0,l), Sig12)
test_stat23 <- mahalanobis(sqrt(N)*(Th[ind2] - Th[ind3]), rep(0,l), Sig23)

# Monte Carlo replicates (for fun)
ts_rep12 <- mahalanobis(mvtnorm::rmvnorm(50000, rep(0, l), Sig12), rep(0,l), Sig12)
ts_rep23 <- mahalanobis(mvtnorm::rmvnorm(50000, rep(0, l), Sig23), rep(0,l), Sig23)

xx <- seq(0,200,.1); yy <- dchisq(xx, l)
par(mfrow=c(2,1), mar=c(2,2,1,1))
hist(ts_rep12, probability=T, breaks=100)
abline(v=test_stat12, col=2, lty=2)
lines(xx, yy, col=4, lwd=2)
hist(ts_rep23, probability=T, breaks=100)
lines(xx, yy, col=4, lwd=2)
abline(v=test_stat23, col=2, lty=2)

# Only case not rejected: Toronto-Ottawa with the first coefficient removed.
pval12 <- c(pchisq(test_stat12, df=l, lower.tail = F), mean(test_stat12 <= ts_rep12))
pval23 <- c(pchisq(test_stat23, df=l, lower.tail = F), mean(test_stat23 <= ts_rep23))
round(rbind(pval12, pval23),4)
#





# Extra extra stuff -------------------------------------------------------
# it takes a bit of time to run...

#########################################################
#### Cautionary tale: pairwise taus between stations ####
#########################################################

# concordance
ijs <- t(combn(1:q,2))
C <- taujack_ms(sapply(Ys, \(Y) Y[-na_rows0,1]), returnC = T)
C <- apply(ijs, 1, \(ij) C[, ij[1], ij[2]]) |> C2Cna(na_rows=na_rows0)

# number of terms of the zeta sum that we want to compute
K <- 1000
list2env(Cna2tsz(C, K, pw=T), envir = environment()) # Th, Sh, and Zh

# quick check
stns_name; cbind(ijs, Th); Sh
matrixcalc::is.positive.definite(Sh)

# select number of terms in the zeta sum
Zhcs <- 4*apply(Zh, c(1,2), \(z) z[1] + cumsum(z[-1])) |> aperm(perm = c(2,3,1))

# plots of the sequences zeta_j (j=1,2,...)
par(mfrow = c(2,3), mar=c(2,2,3,1))
for(i in 1:q) for(j in i:q)
  plot(Zhcs[i,j,-1], main = paste0("pair: (",i,",",j,")"), type="o", cex=.25)

# plots of the sequences zeta_j (j=1,2,...), after cumulative sum is taken
for(i in 1:q) for(j in i:q)
  Zhcs[i,j,] |> plot(main = paste0("pair: (",i,",",j,")"),
                     type="o", cex=.25, ylim = range(Zhcs[1:q,1:q,]))

# test of equality will reject for certain... 
# but it's clear there is a problem with the zetas.
