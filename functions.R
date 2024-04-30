
# functions for the simulation study --------------------------------------
# Knight's extended algorithm
taujack_ms <- function(X, K_serial = 0, returnC = F, flatten = T){

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
  
  if(!flatten){
    if(returnC) return(C)
    if(K_serial == 0) return(list(
      Th = apply(C, c(2,3), sum)/choose(n,2) - 1,
      Sh = 16*cov(apply(ijs, 1, function(k) C[,k[1],k[2]]))/(n*(n-1)),
      ijs = ijs))
    if(K_serial > 0) stop("not implemented")
  }
  
  C <- apply(ijs, 1, function(k) C[,k[1],k[2]])
  if(returnC) return(C)
  
  if(K_serial == 0) return(list(
    Th = apply(C, 2, sum)/choose(n,2) - 1,
    Sh = 16*cov(C)/(n*(n-1)),
    ijs = ijs))
  
  dd <- d-1
  gs <- 2*C/(n-1) - 1
  # Zh <- sapply(0:K_serial, \(k) (n-k-1)/(n-k) * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F]), simplify = "array")
  Zh <- sapply(0:K_serial, \(k) (n-k-1)/n * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F]), simplify = "array")
  if(dd == 1) Zh <- array(Zh, c(1,1,K_serial+1))
  Zh[,,-1] <- Zh[,,-1] + aperm(Zh[,,-1], c(2,1,3))
  Zsum <- apply(Zh[,,-1], c(1,2), sum, na.rm=T)

  return(list(Th = colMeans(gs), Sh = 4*(Zh[,,1] + Zsum), Zh = Zh))
}

# divide-and-conquer algorithm
taujack_dac <- function(X, thresh = 25L, brute_force = F, K_serial = 0, returnC = F){
  
  if(brute_force){
    C <- bruteForce(X, seq = T)
  }else{
    # X <- apply(X, 2, rank)
    C <- dac_seq(X, thresh = thresh, brute_force = F)
  }
  
  if(returnC) return(C)
  
  n <- nrow(X)
  d <- ncol(X)
  dd <- d-1

  gs <- t((2^(1:dd)*t(C)/(n-1) - 1)/(2^(1:dd) - 1))
  Zh <- sapply(0:K_serial, \(k){
    # (n-k-1)/(n-k) * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
    (n-k-1)/n * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
  }, simplify = "array")
  if(dd == 1) Zh <- array(Zh, c(1,1,K_serial+1))
  
  if(K_serial == 0){
    return(list(Th = colMeans(gs),
                Sh = 4*Zh[,,1]))
  }
  
  Zh[,,-1] <- Zh[,,-1] + aperm(Zh[,,-1], c(2,1,3))
  Zsum <- apply(Zh[,,-1], c(1,2), sum, na.rm=T)

  return(list(Th = colMeans(gs),
              Sh = 4*(Zh[,,1] + Zsum),
              Zh = Zh))
}


# divide-and-conquer algorithm
taujack_dac_tune <- function(X, thresh = 25L, brute_force = F, K_serial = 0, returnC = F, tune = 0L){
  
  if(brute_force){
    C <- bruteForce(X, seq = T)
  }else{
    # X <- apply(X, 2, rank)
    C <- dac_seq_tune(X, thresh = thresh, brute_force = F, tune = tune)
  }
  
  if(returnC) return(C)
  
  n <- nrow(X)
  d <- ncol(X)
  dd <- d-1
  
  gs <- t((2^(1:dd)*t(C)/(n-1) - 1)/(2^(1:dd) - 1))
  Zh <- sapply(0:K_serial, \(k){
    # (n-k-1)/(n-k) * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
    (n-k-1)/n * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
  }, simplify = "array")
  if(dd == 1) Zh <- array(Zh, c(1,1,K_serial+1))
  
  if(K_serial == 0){
    return(list(Th = colMeans(gs),
                Sh = 4*Zh[,,1]))
  }
  
  Zh[,,-1] <- Zh[,,-1] + aperm(Zh[,,-1], c(2,1,3))
  Zsum <- apply(Zh[,,-1], c(1,2), sum, na.rm=T)
  
  return(list(Th = colMeans(gs),
              Sh = 4*(Zh[,,1] + Zsum),
              Zh = Zh))
}


# taujack_dacc <- function(X, thresh = 25L, brute_force = F, K_serial = 0, returnC = F){
#   
#   if(brute_force){
#     C <- bruteForce(X, seq = T)
#   }else{
#     # X <- apply(X, 2, rank)
#     C <- dac_seqq(X, thresh = thresh, brute_force = F)
#   }
#   
#   if(returnC) return(C)
#   
#   n <- nrow(X)
#   d <- ncol(X)
#   dd <- d-1
#   
#   gs <- t((2^(1:dd)*t(C)/(n-1) - 1)/(2^(1:dd) - 1))
#   Zh <- sapply(0:K_serial, \(k){
#     # (n-k-1)/(n-k) * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
#     (n-k-1)/n * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
#   }, simplify = "array")
#   if(dd == 1) Zh <- array(Zh, c(1,1,K_serial+1))
#   
#   if(K_serial == 0){
#     return(list(Th = colMeans(gs),
#                 Sh = 4*Zh[,,1]))
#   }
#   
#   Zh[,,-1] <- Zh[,,-1] + aperm(Zh[,,-1], c(2,1,3))
#   Zsum <- apply(Zh[,,-1], c(1,2), sum, na.rm=T)
#   
#   return(list(Th = colMeans(gs),
#               Sh = 4*(Zh[,,1] + Zsum),
#               Zh = Zh))
# }


# taujack_dac_parallel <- function(X, thresh = 25L, brute_force = F, K_serial = 0, returnC = F, mc_cores = 2){
# 
#   n <- nrow(X)
#   d <- ncol(X)
#   dd <- d-1
#   
#   if(brute_force){
#     C <- bruteForce(X, seq = T)
#   }else if(mc_cores == 1){
#       C <- dac_seq(X, thresh = thresh, brute_force = F)
#   }else{
#     l <- mc_cores
#     nn <- round(seq(1,n,length.out=l+1))
#     nn <- matrix(c(nn[1], sapply(nn[-c(1,l+1)], \(x) c(x,x+1)), nn[l+1]),
#                  mc_cores, 2, byrow = T)
#     
#     clus <- parallel::makeCluster(l)
#     parallel::clusterExport(clus, list("seq"), envir = environment())
#     parallel::clusterEvalQ(clus, {
#       library(taujackmin)
#     })
#     
#     C <- parallel::parLapply(cl = clus,
#                              X = apply(nn, 1, \(x) X[seq(x[1],x[2]),], simplify = F),
#                              fun = \(X_sub) dac_seq(X_sub, thresh = thresh, brute_force = F)) |>
#       do.call(what = "rbind")
#     
#     ijs <- combn(l,2)
#     Cs <- parallel::parLapply(cl = clus,
#                               X = apply(ijs, 2, \(x) list(X[seq(nn[x[1],1], nn[x[1],2]),],
#                                                           X[seq(nn[x[2],1], nn[x[2],2]),]), simplify = F),
#                               fun = \(X_sub) dac_seq(X_sub[[1]], X_sub[[2]], thresh = thresh, brute_force = F))
#     
#     for(r in 1:ncol(ijs)){
#       x <- ijs[,r]
#       i1 <- nn[x[1],1]; i2 <- nn[x[1],2]; n1 <- i2 - i1 + 1
#       j1 <- nn[x[2],1]; j2 <- nn[x[2],2]; n2 <- j2 - j1 + 1
#       C[seq(i1, i2),] <- C[seq(i1, i2),] + Cs[[r]][seq(1, n1),]
#       C[seq(j1, j2),] <- C[seq(j1, j2),] + Cs[[r]][n1 + seq(1, n2),]
#     }
#     
#     stopCluster(clus)
#   }
#   
#   if(returnC) return(C)
#   
#   gs <- t((2^(1:dd)*t(C)/(n-1) - 1)/(2^(1:dd) - 1))
#   Zh <- sapply(0:K_serial, \(k){
#     # (n-k-1)/(n-k) * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
#     (n-k-1)/n * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
#   }, simplify = "array")
#   if(dd == 1) Zh <- array(Zh, c(1,1,K_serial+1))
#   
#   if(K_serial == 0){
#     return(list(Th = colMeans(gs),
#                 Sh = 4*Zh[,,1]))
#   }
#   
#   Zh[,,-1] <- Zh[,,-1] + aperm(Zh[,,-1], c(2,1,3))
#   Zsum <- apply(Zh[,,-1], c(1,2), sum, na.rm=T)
#   
#   return(list(Th = colMeans(gs),
#               Sh = 4*(Zh[,,1] + Zsum),
#               Zh = Zh))
# }


# brute force algorithm
taujack_bf <- function(X, seq = F, K_serial = 0, returnC = F){
  C <- bruteForce(X, seq = seq)
  
  if(returnC) return(C)
  
  n <- nrow(X)
  d <- ncol(X)
  if(seq){
    dd <- 1:(d-1)
  }else{
    dd <- d-1
  }
  

  gs <- t((2^dd*t(C)/(n-1) - 1)/(2^dd - 1))
  Zh <- sapply(0:K_serial, \(k){
    # (n-k-1)/(n-k) * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
    (n-k-1)/n * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
  }, simplify = "array")
  if(length(dd) == 1) Zh <- array(Zh, c(1,1,K_serial+1))
  
  if(K_serial == 0){
    return(list(Th = colMeans(gs),
                Sh = 4*Zh[,,1]))
  }
  
  Zh[,,-1] <- Zh[,,-1] + aperm(Zh[,,-1], c(2,1,3))
  Zsum <- apply(Zh[,,-1], c(1,2), sum, na.rm=T)

  return(list(Th = colMeans(gs),
              Sh = 4*(Zh[,,1] + Zsum),
              Zh = Zh))
}


taujack_bf_parallel <- function(X, seq = F, K_serial = 0, returnC = F, mc_cores = 2){
  n <- nrow(X)
  d <- ncol(X)
  dd <- d-1
  
  if(mc_cores == 1){
    C <- bruteForce(X, seq = seq)
  }else{
    l <- mc_cores
    nn <- round(seq(1,n,length.out=l+1))
    nn <- matrix(c(nn[1], sapply(nn[-c(1,l+1)], \(x) c(x,x+1)), nn[l+1]),
                 mc_cores, 2, byrow = T)
    
    clus <- parallel::makeCluster(l)
    parallel::clusterExport(clus, list("seq"), envir = environment())
    parallel::clusterEvalQ(clus, {
      library(taujackmin)
    })
    
    C <- parallel::parLapply(cl = clus,
                             X = apply(nn, 1, \(x) X[seq(x[1],x[2]),], simplify = F),
                             fun = \(X_sub) bruteForce(X_sub, seq)) |>
      do.call(what = "rbind")
    
    ijs <- combn(l,2)
    Cs <- parallel::parLapply(cl = clus,
                              X = apply(ijs, 2, \(x) list(X[seq(nn[x[1],1], nn[x[1],2]),],
                                                          X[seq(nn[x[2],1], nn[x[2],2]),]), simplify = F),
                              fun = \(X_sub) bruteForce2(X_sub[[1]], X_sub[[2]], seq))
    
    for(r in 1:ncol(ijs)){
      x <- ijs[,r]
      i1 <- nn[x[1],1]; i2 <- nn[x[1],2]; n1 <- i2 - i1 + 1
      j1 <- nn[x[2],1]; j2 <- nn[x[2],2]; n2 <- j2 - j1 + 1
      C[seq(i1, i2),] <- C[seq(i1, i2),] + Cs[[r]][seq(1, n1),]
      C[seq(j1, j2),] <- C[seq(j1, j2),] + Cs[[r]][n1 + seq(1, n2),]
    }
    
    stopCluster(clus)
  }
  
  if(returnC) return(C)
  
  gs <- t((2^(1:dd)*t(C)/(n-1) - 1)/(2^(1:dd) - 1))
  Zh <- sapply(0:K_serial, \(k){
    # (n-k-1)/(n-k) * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
    (n-k-1)/n * cov(gs[1:(n-k),,drop=F], gs[(k+1):n,,drop=F])
  }, simplify = "array")
  if(dd == 1) Zh <- array(Zh, c(1,1,K_serial+1))
  
  if(K_serial == 0){
    return(list(Th = colMeans(gs),
                Sh = 4*Zh[,,1]))
  }
  
  Zh[,,-1] <- Zh[,,-1] + aperm(Zh[,,-1], c(2,1,3))
  Zsum <- apply(Zh[,,-1], c(1,2), sum, na.rm=T)
  
  return(list(Th = colMeans(gs),
              Sh = 4*(Zh[,,1] + Zsum),
              Zh = Zh))
}



# Helpers for application -------------------------------------------------
C2Cna <- function(C, na_rows){
  C_temp <- matrix(NA, nrow(C) + length(na_rows), ncol(C))
  C_temp[-na_rows,] <- C
  C_temp
}

Cna2tsz <- function(Cna, K_serial = 0, pw = F){
  cat("performing Cna2tsz:")
  N <- nrow(Cna)
  n <- sum(apply(!is.na(Cna),1,all))
  dd <- ifelse(pw, 1, ncol(Cna))

  gh <- t((2^(1:dd)*t(Cna)/(n-1) - 1)/(2^(1:dd) - 1))

  Zh <- sapply(0:K_serial, \(k){
    if(k %% 5 == 0) cat(round(100*k/(K_serial + 1),1), "% -- ")
    
    G1 <- gh[1:(N-k),,drop=F]
    G2 <- gh[(k+1):N,,drop=F]
    nn <- sum(apply(!is.na(G1*G2),1,all)) # let us use only "complete" obs.
    
    # (nn/n) * (nn-1)/nn * cov(G1, G2, use = "na.or.complete")
    (nn-1)/n * cov(G1, G2, use = "na.or.complete")
  }, simplify = "array")
  
  if(length(dim(Zh)) < 3) Zh <- array(Zh, c(1,1,K_serial+1))
  if(K_serial > 0){
    Zh[,,-1] <- Zh[,,-1] + aperm(Zh[,,-1], c(2,1,3))
    Zsum <- apply(Zh[,,-1], c(1,2), sum, na.rm=T)
  }else{
    Zsum <- 0
  }
  
  cat("\n")
  return(list(Th = colMeans(gh, na.rm=T),
              Sh = 4*(Zh[,,1] + Zsum),
              Zh = Zh,
              gh = gh))
}

Zh2Sh <- function(Zh, K=0){
  if(K > 0){
    Zsum <- apply(Zh[,,1 + 1:K], c(1,2), sum)
  }else{
    Zsum <- 0
  }
  return(Sh = 4*(Zh[,,1] + Zsum))
}
