# Informal tests for checking the recurence relations for tau (odd dimension)


# packages ----------------------------------------------------------------
library(tautests) # (not available online)




# concordance computations ------------------------------------------------

# dimension 3
d <- 3
n <- 10
X <- replicate(d, sample(n))

C0 <- sapply(1:n, function(i){sum(apply(X[i,] < t(X[-i,]), 2, all)|apply(X[i,] > t(X[-i,]), 2, all))})
# C0; conc(X)

C2 <- tautests::conc(X)
C3 <- n-1- ((n - 1 - C2[,1,2]) + (n - 1 - C2[,2,3]) + (n - 1 - C2[,1,3]))/2
C0; C3

D2 <- n-1-C2
D3 <- n-1-C3
((D2[,1,2]) + (D2[,2,3]) + (D2[,1,3]))/2

T2 <- -1+2*C2/(n-1)
T3 <- (-1+4*C3/(n-1))/3

rowSums(cbind(T2[,1,2],T2[,2,3],T2[,1,3]))/3
T3

sum(C0)/2
sum(C3)/2



# dimension 5
d <- 5
n <- 50
X <- replicate(d, sample(n))

C0 <- sapply(1:n, function(i){sum(apply(X[i,] < t(X[-i,]), 2, all)|apply(X[i,] > t(X[-i,]), 2, all))})
C0
T0 <- (-1+2^(d-1)*C0/(n-1))/(2^(d-1)-1)
# C0; conc(X)

C2 <- tautests::conc(X)
com <- combn(1:d, 3)

C3 <- matrix(nrow=n, ncol=0)
for(k in 1:ncol(com)){
  C3 <- cbind(C3,
              n-1- ((n - 1 - C2[,com[1,k],com[2,k]]) + (n - 1 - C2[,com[2,k],com[3,k]]) + (n - 1 - C2[,com[1,k],com[3,k]]))/2
  )
}
C2 <- apply(combn(d,2), 2, \(ij) C2[,ij[1], ij[2]])

C4 <- apply(combn(d, 4), 2, function(com){
  sapply(1:n, function(i){sum(apply(X[i,com] < t(X[-i,com]), 2, all)|apply(X[i,com] > t(X[-i,com]), 2, all))})
})

T2 <- -1+2*C2/(n-1)
T3 <- (-1+4*C3/(n-1))/3
T4 <- (-1+8*C4/(n-1))/7
T5 <- 2^(d-1)/(2^(d-1)-1) * ((-1/2)^2 * rowSums(T2) +
                               (-1/2)^3 * (2^2-1) * rowSums(T3) +
                               (-1/2)^4 * (2^3-1) * rowSums(T4))

rowSums(T2[,c(1,2,5)])/3 - T3[,1]
T0 - T5

