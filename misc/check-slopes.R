# Script that generates the plots for the benchmark analysis

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


times0[fun == "BF" & n %in% 2^(14:16) & type == "median runtime", .(n, p, time, tau)] |>
  dcast(n+tau~p, value.var = "time")


times0[order(fun, p, n)][type == "median runtime" & tau == .5 & fun == "DAC", 
                         median(log(time,2)), .(fun, p, n)]

tt <- times0[order(fun, p, n)][type == "median runtime" & tau == .25 & fun == "DAC", 
                         median(log(time,2)), .(p, n)] |>
  dcast(n~p)

apply(tt, 2, diff)

t(apply(apply(tt, 2, diff)[,2:5],1,diff))
