times0[order(fun, p, n)][type == "median runtime" & tau == .5 & fun == "DAC", 
                         median(log(time,2)), .(fun, p, n)]

tt <- times0[order(fun, p, n)][type == "median runtime" & tau == .5 & fun == "DAC", 
                         median(log(time,2)), .(p, n)] |>
  dcast(n~p)

apply(tt, 2, diff)
