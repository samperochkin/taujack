
Pm <- NULL
Ps <- NULL
Pn <- NULL
Pt <- NULL

for(h in 0:23){
  P <- readRDS(paste0("app/P_", h, ".rds"))
  
  Pm <- rbind(Pm, P[,1])
  Ps <- rbind(Ps, P[,2])
  Pn <- rbind(Pn, P[,3])
  Pt <- rbind(Pt, P[,4])
}

B <- lapply(1:10, \(i) cbind(sin(2*i*pi*(0:23)/24), cos(2*i*pi*(0:23)/24))) |> 
  do.call(what="cbind") |> cbind(1)

Pm <- B %*% (MASS::ginv(B) %*% Pm)
Ps <- B %*% (MASS::ginv(B) %*% Ps)
Pn <- B %*% (MASS::ginv(B) %*% Pn)
Pt <- B %*% (MASS::ginv(B) %*% Pt)

data0$res <- as.numeric(NA)
for(h in 0:23){
  indh <- ind[data0[ind,]$hour == h]
  X <- model.matrix(ff, data0[indh,])
  mu <- X %*% Pm[h+1,]
  sig <- exp(X %*% Ps[h+1,])
  nu <- X %*% Pn[h+1,]
  tau <- exp(X %*% Pt[h+1,])
  
  res <- gamlss.dist::pST2(data0[indh,]$temp,mu,sig,nu,tau)
  hist(res, probability = T, breaks=20)
  data0[indh,]$res <- res
}

hist(data0$res, probability = T, breaks=20)

# quick check on about 9 years of data (so that the series is uninterupted)
# conclusion: small correlation over few hours, but quickly vanishes.
# this is satisfactory
acf(data0$res[25:85000], lag.max = 48)
acf(data0$res[25:85000], lag.max = 48, ylim = c(-.1,.1))
acf(data0$res[25:85000], lag.max = 24*30, ylim = c(-.025,.025))
acf(data0$res[25:85000], lag.max = 24*30*12*2, ylim = c(-.025,.025))

pacf(data0$res[25:85000], lag.max = 48)
pacf(data0$res[25:85000], lag.max = 48, ylim = c(-.1,.1))
pacf(data0$res[25:85000], lag.max = 24*30, ylim = c(-.025,.025))
pacf(data0$res[25:85000], lag.max = 24*30*12*2, ylim = c(-.025,.025))


nrow(data0)
which(is.na(data0$res))
