# Script that generates the plots for the benchmark analysis

# packages ----------------------------------------------------------------
library(ggplot2)
library(data.table)



# load results and format -------------------------------------------------
times <- list.files("benchmark-low-dim/", pattern = "times[1-9]", full.names = TRUE) |> lapply(fread) |> rbindlist()
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
type0 <- "median runtime"
tau0 <- 1
fun0 <- "DAC"
p0 <- 2
times00 <- times0[order(fun, n)][tau == tau0 & type == type0 & p==p0 & fun == fun0]
diff(times00$ltime)

p0 <- 2
times00 <- times0[order(fun, n)][tau == tau0 & type == type0 & p==p0 & fun == fun0]
diff(times00$ltime)

fun0 <- "KO"
p0 <- 2
times00 <- times0[order(fun, n)][tau == tau0 & type == type0 & p==p0 & fun == fun0]
diff(times00$ltime)

fun0 <- "KE"
p0 <- 2
times00 <- times0[order(fun, n)][tau == tau0 & type == type0 & p==p0 & fun == fun0]
diff(times00$ltime)

# plots -------------------------------------------------------------------

# for main document
gg <- ggplot(times0[p %in% c(2,4,6,10) & type != "mean runtime"],
             aes(x = log(n,2), y = log(time,2), linetype=as.factor(p), col=as.factor(fun))) +
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
  scale_x_continuous(breaks = seq(6,34,2)) +
  scale_y_continuous(breaks = seq(-14,22,4)) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](average~runtime))) +
  labs(col="algorithm") +
  scale_linetype_manual(name = "dimension (p)", values = c(1,2,3,4)) +
  geom_line(size=.5) +
  geom_point(size=2) +
  facet_grid(tau~type, labeller = label_bquote(rows = "t" == .(tau)))
gg



# for appendix
gg <- ggplot(times0[p %in% c(2,4,6,10) & type == "median runtime" & tau %in% c(0,.5,1)],
             aes(x = log(n,2), y = log(time,2), linetype=as.factor(p), col=as.factor(fun))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=9),
        legend.position = "right",
        legend.direction = "vertical",
        legend.background = element_rect(colour="black"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size=10, margin = margin(r = 10, unit = "pt")),
        # legend.spacing.x = unit(.5, "cm"),
        panel.spacing = unit(0.5, "lines"),
        strip.text.x = element_text(size = 10)) +
  guides(col = guide_legend(ncol=2), linetype = guide_legend(ncol=2)) +
  scale_x_continuous(breaks = seq(6,30,2)) +
  scale_y_continuous(breaks = seq(-14,10,4)) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](median~runtime))) +
  labs(col="algorithm") +
  scale_linetype_manual(name = "dimension (p)", values = c(1,2,3,4)) +
  geom_line(size=.5) +
  facet_grid(~tau, labeller = label_bquote(cols = "t" == .(tau)))
gg


gg <- ggplot(times0[p %in% c(2,4,6,10) & type == "median runtime" & fun %in% c("KO", "KE", "DAC") & tau == 1],
             aes(x = n, y = time, linetype=as.factor(p), col=as.factor(fun))) +
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
  scale_y_continuous(breaks = seq(-14,6,4)) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](average~runtime))) +
  labs(col="algorithm") +
  scale_linetype_manual(name = "dimension (p)", values = c(1,2,3,4)) +
  geom_line(size=.5) +
  geom_point(size=2) +
  facet_grid(tau~type, labeller = label_bquote(rows = "t" == .(tau)))
gg



# stuff to compare with
nn <- seq(1,1e3, length.out=10)
k <- 4
Y <- sapply(1:k, \(x) x*nn*log(nn,2))
Z <- sapply(1:k, \(x) nn*log(nn,2)^x)

plot(nn, Y[,1], type="o", ylim=range(Y), pch=19)
for(i in 2:k) lines(nn, Y[,i], col = i, type="o", pch=19)

plot(nn, Z[,1], type="o", ylim=range(Z), pch=19)
for(i in 2:k) lines(nn, Z[,i], col = i, type="o", pch=19)

plot(nn, Z[,4]-Z[,2], type="o", ylim=range(Z), pch=19)
lines(nn, Y[,4]-Y[,2], type="o", pch=19)

# this might be the best thing
plot(nn, Z[,4]-Z[,2], type="o", pch=19)
plot(nn, Y[,4]-Y[,2], type="o", pch=19)



# it seems that the f(p)*nn*log(nn,2) is the more plausible
times1 <- times0[p %in% c(2,4) & type == "median runtime" & fun %in% c("DAC") & tau == 1]
times1 <- dcast(times1, n~p, value.var = "time")
names(times1) <- c("n", "t2", "t4")
times1[, diff := t4-t2]

gg <- ggplot(times1,
             aes(x = n, y = diff)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=9),
        strip.text.x = element_text(size = 10)) +
  geom_line(size=.5) +
  geom_point(size=2)
gg
