# Script that generates the plots for the benchmark analysis

# packages ----------------------------------------------------------------
library(ggplot2)
library(data.table)



# load results and format -------------------------------------------------
times <- c(
    list.files("benchmark-p2", pattern = "times_0_29_[1-9]", full.names = TRUE),
    list.files("benchmark-p4-0", pattern = "times[1-9]", full.names = TRUE),
    list.files("benchmark-p4", pattern = "times[1-9]", full.names = TRUE)
  ) |> lapply(fread) |> rbindlist()
times <- times[order(p,tau,n)]

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




# Quick check of slopes ---------------------------------------------------

times0[, ltime := log(time,2)]
times0[, slope := c(diff(ltime),NA), .(p,tau,type)]

# plots -------------------------------------------------------------------

gg <- ggplot(times0[type == "median runtime"],
             aes(x = log(n,2), y = log(time,2), col=slope)) +
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
  scale_x_continuous(breaks = seq(6,34,4)) +
  scale_y_continuous(breaks = seq(-14,22,4)) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](median~runtime))) +
  viridis::scale_color_viridis(name = "slope", breaks = c(.2,.6,1,1.4)) +
  # viridis::scale_color_viridis(option = "H", name = "slope", breaks = c(.2,.6,1,1.4)) +
  geom_line(size=.75) +
  geom_point(size=1.25) +
  facet_wrap(~p+tau, labeller = label_bquote(rows = "t = "*.(tau)*",  p = "*.(p)))
gg

times0[n == 2^26 & p == 2 & type == "median runtime"]
times0[n == 2^24 & p == 4 & tau == 0 & type == "median runtime"]
times0[n == 2^24 & p == 4 & tau == .5 & type == "median runtime"]

pdf("figures/extra.pdf", height = 6.5/2.5, width = 6.5)
print(gg)
dev.off()





