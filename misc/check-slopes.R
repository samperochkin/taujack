# packages ----------------------------------------------------------------
library(ggplot2)
library(data.table)



# load results and format -------------------------------------------------
times <- list.files("benchmark/", pattern = "times[1-9]", full.names = TRUE) |> lapply(fread) |> rbindlist()
times0 <- times[, .(mean_time = mean(time), min_time = min(time), max_time = max(time), N = .N),
                .(n,tau,p,fun)]

times0 <- rbind(
  times[, .(time = mean(time), type = "mean", N = .N), .(n,tau,p,fun)],
  times[, .(time = median(time), type = "median", N = .N), .(n,tau,p,fun)],
  times[, .(time = min(time), type = "min", N = .N), .(n,tau,p,fun)],
  times[, .(time = max(time), type = "max", N = .N), .(n,tau,p,fun)]
)
times0$type <- factor(times0$type,
                      levels = c("min", "mean", "median", "max"),
                      labels = c("min. runtime", "mean runtime", "median runtime", "max. runtime"))
times0[grepl("BF", fun), fun := "BF"]
times0[grepl("KO", fun), fun := "KO"]
times0[grepl("KE", fun), fun := "KE"]
times0[grepl("DAC", fun), fun := "DAC"]



# tau0 <- .5
tau0 <- .75
# times0[order(fun, p, n)][type == "median runtime" & tau == tau0 & fun == "DAC", 
#                          median(log(time,2)), .(fun, p, n)]
tt <- times0[order(fun, p, n)][type == "median runtime" & tau == tau0 & fun == "DAC", 
                         median(log(time,2)), .(p, n)] |> dcast(n~p)
tt <- cbind(tt$n[-1], log(tt$n,2)[-1], apply(tt[,-1], 2, diff))
tt

t(apply(tt[,-(1:2)], 1, diff))
