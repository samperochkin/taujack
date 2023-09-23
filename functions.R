
# functions for the simulation study --------------------------------------
# Knight's extended algorithm
taujack_ms <- function(X, returnC = F){

  n <- nrow(X); d <- ncol(X)
  ijs <- t(sapply(1:choose(d,2), function(r){
    j <- ceiling((1 + sqrt(1 + 8*r))/2)
    i <- r - choose(j-1,2)
    return(c(i,j))
  }))
  
  C <- array(n-1, dim = c(n,d,d))
  
  for (i in 1:(d-1)){
    ord <- order(X[, i])
    for (j in (i+1):d){
      C[ord, i, j] <- C[ord, j, i] <- n - 1 - countSwaps(X[ord,j])
    }
  }
  
  if(returnC) return(C)
  
  return(list(tau = apply(C, c(2,3), sum)/choose(n,2) - 1,
              sigma2 = 16*cov(apply(ijs, 1, function(k) C[,k[1],k[2]]))/(n*(n-1))))
}

# divide-and-conquer algorithm
taujack_dac <- function(X, thresh = 25L, brute_force = F, K_serial = 0, returnC = F){
  
  if(brute_force){
    C <- bruteForce(X, seq = T)
  }else{
    X <- apply(X, 2, rank)
    C <- dac_seq(X, thresh = thresh, brute_force = F)
  }
  
  if(returnC) return(n)
  
  n <- nrow(X)
  d <- ncol(X)
  dd <- d-1

  gs <- t((2^(1:dd)*t(C)/(n-1) - 1)/(2^(1:dd) - 1))
  Zs <- sapply(0:K_serial, \(k) cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F]), simplify = "array")
  if(dd == 1) Zs <- array(Zs, c(1,1,K_serial+1))
  if(K_serial > 0){
    Zsum <- apply(Zs[,,-1], c(1,2), sum)
    Zsum <- Zsum + t(Zsum)
  }else{
    Zsum <- 0
  }

  return(list(tau = colMeans(gs),
              sigma2 = 4*(Zs[,,1] + Zsum),
              zetas = Zs))
}

# brute force algorithm
taujack_bf <- function(X, seq = F, K_serial = 0, returnC = F){
  C <- bruteForce(X, seq = seq)
  
  if(returnC) return(C)
  
  n <- nrow(X)
  d <- ncol(X)
  dd <- d-1

  gs <- t((2^(1:dd)*t(C)/(n-1) - 1)/(2^(1:dd) - 1))
  Zs <- sapply(0:K_serial, \(k) cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F]), simplify = "array")
  if(dd == 1) Zs <- array(Zs, c(1,1,K_serial+1))
  if(K_serial > 0){
    Zsum <- apply(Zs[,,-1], c(1,2), sum)
    Zsum <- Zsum + t(Zsum)
  }else{
    Zsum <- 0
  }
  
  return(list(tau = colMeans(gs),
              sigma2 = 4*(Zs[,,1] + Zsum),
              zetas = Zs))
  
}
