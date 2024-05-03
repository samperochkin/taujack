# Script that generates the plots for the benchmark analysis

# packages ----------------------------------------------------------------
library(ggplot2)
library(data.table)



# load results and format -------------------------------------------------
pat <- "times[1-9]"

times <- list.files("benchmark-p4-unif", pattern = pat, full.names = TRUE) |> lapply(fread) |> rbindlist()
times <- times[order(n)]
times0 <- times[, .(mean_time = mean(time), min_time = min(time), max_time = max(time), N = .N),
                .(n,p,fun)]
times0 <- rbind(
  times[, .(time = mean(time), type = "mean", N = .N), .(n,p,fun)],
  times[, .(time = median(time), type = "median", N = .N), .(n,p,fun)],
  times[, .(time = min(time), type = "min", N = .N), .(n,p,fun)],
  times[, .(time = max(time), type = "max", N = .N), .(n,p,fun)]
)
times0$type <- factor(times0$type,
                      levels = c("min", "mean", "median", "max"),
                      labels = c("min. runtime", "mean runtime", "median runtime", "max. runtime"))
times0[grepl("BF", fun), fun := "BF"]
times0[grepl("KO", fun), fun := "KO"]
times0[grepl("KE", fun), fun := "KE"]
times0[grepl("DAC", fun), fun := "DAC"]




# Quick check of slopes ---------------------------------------------------

times0[, ltime := log(time,2)]
type0 <- "median runtime"
fun0 <- "DAC"
p0 <- 4
times00 <- times0[order(fun, n)][type == type0 & p==p0 & fun == fun0]
round(rbind(log(times00$n[-1],2), diff(times00$ltime)/diff(log(times00$n,2))),2)
plot(times00$time[-nrow(times00)], diff(times00$time)/diff(times00$n))
times00$slope <- c(diff(times00$ltime), NA)/c(diff(log(times00$n,2)), NA)
plot(times00$slope)


ii <- c(round(nrow(times00)/2),nrow(times00))
diff(times00$ltime[ii])/diff(log(times00$n[ii],2))


ii <- c(which.min(abs(times00$n - 2^17)), which.min(abs(times00$n - 2^19)))
diff(times00$ltime[ii])/diff(log(times00$n[ii],2))
1 + p0 * (log(19,2) - log(18,2))

ii <- c(which.min(abs(times00$n - 2^18)), nrow(times00))
diff(times00$ltime[ii])/diff(log(times00$n[ii],2))
1 + p0 * (log(19,2) - log(18,2))


# plots -------------------------------------------------------------------

gg <- ggplot(times00,
             aes(x = log(n,2), y = log(time,2), linetype=as.factor(p), col=slope)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=9),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_rect(colour="black"),
        legend.title = element_text(size = 10),
        legend.text=element_text(size=10),
        panel.spacing = unit(0.5, "lines"),
        strip.text.x = element_text(size = 10)) +
  # guides(col = guide_legend(ncol=2), linetype = guide_legend(ncol=2)) +
  scale_x_continuous(breaks = seq(6,30,2)) +
  scale_y_continuous(breaks = seq(-14,22,2)) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](median~runtime))) +
  labs(col="algorithm") +
  scale_linetype_manual(name = "dimension (p)", values = c(1,2,3,4)) +
  viridis::scale_color_viridis(name = "slope") +
  geom_line(size=.5) +
  geom_point(size=2) +
  facet_wrap(~type)
gg


pdf("figures/extra_p4_unif.pdf", height = 6.5/3, width = 6.5/3)
print(gg)
dev.off()
