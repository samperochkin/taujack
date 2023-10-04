library(ggplot2)
library(data.table)

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


type0 <- "max. runtime"
type0 <- "median runtime"
times0[order(fun, n, p)][type == type0 & tau == 1, diff(log(time,2)), .(fun, p)]
times0[order(fun, n, p)][type == type0 & tau == 1 & n %in% 2^(16:21), diff(log(time,2)), .(fun, p)]
times0[order(fun, n, p)][type == type0 & tau == 1 & n %in% 2^(16:18), diff(log(time,2)), .(fun, p)]
times0[order(fun, n, p)][type == type0 & tau == 1 & n %in% 2^(19:21), diff(log(time,2)), .(fun, p)]
times0[order(fun, n, p)][type == type0 & tau == 1 & n %in% 2^(17:21), diff(log(time,2)), .(fun, p)]
times0[order(fun, p, n)][type == "median runtime" & tau == 1, median(log(time,2)), .(fun, p, n)]


tt <- times0[order(fun, n, p)][type == type0 & tau == 1, diff(log(time,2)), .(fun, p)]
plot(tt[fun == "DAC" & p == 2]$V1, type = "o")
lines(tt[fun == "DAC" & p == 4]$V1, type = "o", col=2)
lines(tt[fun == "DAC" & p == 6]$V1, type = "o", col=3)
lines(tt[fun == "DAC" & p == 10]$V1, type = "o", col=3)


gg <- ggplot(times0[p %in% c(2,4,6,10) & type != "mean runtime"],
             aes(x = log(n,2), y = log(time,2), linetype=as.factor(p), col=fun)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=9),
        # legend.position = "bottom",
        # legend.direction = "horizontal",
        # legend.box = "vertical",
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.background = element_rect(colour="gray"),
        legend.title = element_text(size = 10),
        legend.text=element_text(size=10),
        panel.spacing = unit(0.5, "lines"),
        strip.text.x = element_text(size = 10)) +
  # guides(col = guide_legend(ncol=2), linetype = guide_legend(ncol=2)) +
  scale_y_continuous(breaks = seq(-14,6,4)) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](average~runtime))) +
  labs(col="algorithm") +
  scale_linetype_manual(name = "dimension (p)", values = c(1,2,3,4)) +
  geom_line(size=.5) +
  # geom_point(size=.25) +
  # geom_line() +
  # geom_point() +
  # facet_grid(tau~type, labeller = label_bquote(cols = tau == .(tau)))
  # facet_grid(tau~type)
  # facet_grid(type~tau, labeller = label_bquote(cols = tau == .(tau)))
  facet_grid(tau~type, labeller = label_bquote(rows = tau == .(tau)))
gg

gg <- ggplot(times0[p %in% c(2,4,6,10) & type == "median runtime" & tau %in% c(0,.5,1)],
             aes(x = log(n,2), y = log(time,2), linetype=as.factor(p), col=fun)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=9),
        # legend.position = "bottom",
        # legend.box = "horizontal",
        # legend.direction = "vertical",
        legend.position = "right",
        legend.direction = "vertical",
        legend.background = element_rect(colour="gray"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10, margin = margin(r = 10, unit = "pt")),
        # legend.spacing.x = unit(.5, "cm"),
        panel.spacing = unit(0.5, "lines"),
        strip.text.x = element_text(size = 10)) +
  guides(col = guide_legend(ncol=2), linetype = guide_legend(ncol=2)) +
  # guides(linetype = guide_legend(ncol=2)) +
  # guides(colour = "none") +
  scale_x_continuous(breaks = seq(6,20,2)) +
  scale_y_continuous(breaks = seq(-14,2,4)) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](median~runtime))) +
  labs(col="algorithm") +
  scale_linetype_manual(name = "dimension (p)", values = c(1,2,3,4)) +
  geom_line(size=.5) +
  # geom_point(size=.25) +
  facet_grid(~tau, labeller = label_bquote(cols = tau == .(tau)))
gg

