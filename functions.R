
# Knight's extended algorithm
taujack_ms <- function(X){

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
  
  return(list(tau = apply(C, c(2,3), sum)/choose(n,2) - 1,
              sigma2 = 16*cov(apply(ijs, 1, function(k) C[,k[1],k[2]]))/(n*(n-1))))
}

# divide-and-conquer algorithm
taujack_dac <- function(X, thresh = 100L, brute_force = F){
  
  n <- nrow(X)
  if(brute_force){
    C <- bruteForce(X)
  }else{
    X <- apply(X, 2, rank)
    C <- dac(X, thresh = thresh, brute_force = brute_force)
  }

  return(list(tau = sum(C)/choose(n,2) - 1,
              sigma2 = 16*var(C)/(n*(n-1))))
}

# brute force algorithm
taujack_bf <- function(X){
  
  n <- nrow(X)
  C <- bruteForce(X)
  return(list(tau = sum(C)/choose(n,2) - 1,
              sigma2 = 16*var(C)/(n*(n-1))))
}
