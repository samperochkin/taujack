# Checking formulas for the number of pairs of groups with certain 
# relations (e.g., concordant, discordant, concordant on all but k 
# dimensions and unknown for the others, ...)
# These are not given in the paper.


# first investigation (simple) --------------------------------------------
countCases <- function(d){
  V <- expand.grid(lapply(1:d, \(x) 0:1))
  cases <- sapply(1:nrow(V), \(k){
    print(k)
    sapply(k:nrow(V), \(l){
      v <- V[k,] - V[l,]
      if(all(v == 0)){
        a <- c(1,0,0,0,0)
      }else if(all(range(v) == c(-1,1))){
        a <- c(0,0,1,0,0)
      }else if(all(v == 1) | all(v == -1)){
        a <- c(0,1,0,0,0)
      }else if(sum(v != 0) == 1){
        a <- c(0,0,0,1,0)
      }else{
        a <- c(0,0,0,0,1)
      }
      a
    })
  }) |> do.call(what = "cbind") |> t()
  
  colSums(cases)
}

xx <- 2:8
cases <- sapply(xx, countCases)
cases
colSums(cases)

cases
rbind(2^xx, 1, 2^(xx-1)*(2^xx + 1) - 3^xx, 2^(xx-1)*xx, 3^xx - 2^xx - 1 - 2^(xx-1)*xx)
cases - rbind(2^xx, 1, 2^(xx-1)*(2^xx + 1) - 3^xx, 2^(xx-1)*xx, 3^xx - 2^xx - 1 - 2^(xx-1)*xx)

xx <- 3:50
crits <- cbind(log(2*(3^(xx) - 2^(xx - 1) - 1), 2^(xx)),
               log(2*(3^(xx) - (xx + 1)*2^(xx - 1) - 1), 2^(xx)))

plot(xx, crits[,1], type= "o", ylim=range(crits), pch=19, cex=.5)
lines(xx, crits[,2], type= "o", col=3, pch=19, cex=.5)



# More cases --------------------------------------------------------------
d <- 4
V <- expand.grid(lapply(1:d, \(x) 0:1))
cases <- sapply(1:nrow(V), \(k){
  print(k)
  sapply(k:nrow(V), \(l){
    v <- V[k,] - V[l,]
    if(all(v == 0)){
      a <- c(1,0,0,0,0,0,0)
    }else if(all(range(v) == c(-1,1))){
      a <- c(0,0,1,0,0,0,0)
    }else if(all(v == 1) | all(v == -1)){
      a <- c(0,1,0,0,0,0,0)
    }else if(sum(v != 0) == 1){
      a <- c(0,0,0,1,0,0,0)
    }else if(sum(v != 0) == 2){
      a <- c(0,0,0,0,1,0,0)
    }else if(sum(v != 0) == 3){
      a <- c(0,0,0,0,0,1,0)
    }else{
      # overflow cases
      a <- c(0,0,0,0,0,0,1)
    }
    a
  })
}) |> do.call(what = "cbind") |> t()


cases
colSums(cases)
c(2^d, 1, 2^(d-1)*(2^d + 1) - 3^d, 2^(d-1)*d, 2^(d-3)*(d-1)*d, 2^(d-4)*(d-2)*(d-1)*d/3)
colSums(cases)[1:6] - c(2^d, 1, 2^(d-1)*(2^d + 1) - 3^d, 2^(d-1)*d, 2^(d-3)*(d-1)*d, 2^(d-4)*(d-2)*(d-1)*d/3)
