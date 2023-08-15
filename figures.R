library(ggplot2)
library(data.table)

times <- list.files("benchmark/", pattern = "times[1-9]", full.names = TRUE) |> lapply(fread) |> rbindlist()
times0 <- times[, .(mean_time = mean(time), min_time = min(time), max_time = max(time), N = .N),
                  .(n,tau,p,fun)]
# times0 <- times0[!(fun == "d-a-c" & p %in% c(4,8))]
# times0 <- times0[!(fun == "d-a-c" & p != 8)]

k <- length(unique(times0$p))

gg <- ggplot(times0[p %in% c(2,4,6,10)], aes(x = log(n,2), y = log(mean_time,2), linetype=as.factor(p), col=as.factor(fun))) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=11),
        legend.background = element_rect(colour="gray"),
        legend.title = element_text(size = 10),
        legend.text=element_text(size=10),
        panel.spacing = unit(0.5, "lines"),
        strip.text.x = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(-14,0,2)) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](average~runtime))) +
  labs(col="algorithm") +
  scale_linetype_manual(name = "dimension (p)", values = c(1,2,3,4)) +
  geom_line(linewidth=.5) +
  geom_point(size=.25) +
  geom_line() +
  geom_point() +
  facet_wrap(~tau, labeller = label_bquote(cols = tau == .(tau)))
gg


pdf("res.pdf", height = 3.5, width = 3*3.5)
print(gg)
dev.off()




gg <- ggplot(times0, aes(x = log(n,2), y = log(min_time,2), linetype=as.factor(p), col=as.factor(fun))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text=element_text(size=11),
        legend.background = element_rect(colour="gray"),
        legend.title = element_text(size = 10),
        legend.text=element_text(size=10),
        panel.spacing = unit(0.5, "lines"),
        strip.text.x = element_text(size = 10)) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](minimum~runtime))) +
  labs(col="algorithm") +
  scale_linetype_manual(name = "dimension (p)", values = 1:k) +
  geom_line(size=.25) +
  geom_point(size=.25) +
  geom_line() +
  geom_point() +
  facet_wrap(~tau, labeller = label_bquote(cols = tau == .(tau)))
gg



gg <- ggplot(times0, aes(x = log(n,2), y = log(max_time,2), linetype=as.factor(p), col=as.factor(fun))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text=element_text(size=11),
        legend.background = element_rect(colour="gray"),
        legend.title = element_text(size = 10),
        legend.text=element_text(size=10),
        panel.spacing = unit(0.5, "lines"),
        strip.text.x = element_text(size = 10)) +
  # coord_fixed(ratio = 1) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](maximum~runtime))) +
  labs(col="algorithm") +
  scale_linetype_manual(name = "dimension (p)", values = 1:k) +
  geom_line(size=.25) +
  geom_point(size=.25) +
  geom_line() +
  geom_point() +
  facet_wrap(~tau, labeller = label_bquote(cols = tau == .(tau)))
gg




gg <- ggplot(times0, aes(x = n, y = mean_time, linetype=as.factor(p), col=as.factor(fun))) +
  # gg <- ggplot(times0[fun != "brute force"], aes(x = n, y = mean_time, linetype=as.factor(p), col=as.factor(fun))) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text=element_text(size=11),
        legend.background = element_rect(colour="gray"),
        legend.title = element_text(size = 10),
        legend.text=element_text(size=10),
        panel.spacing = unit(0.5, "lines"),
        strip.text.x = element_text(size = 10)) +
  xlab(expression(log[2]~n)) +
  ylab(expression(log[2](average~runtime))) +
  labs(col="algorithm") +
  scale_linetype_manual(name = "dimension (p)", values = 1:k) +
  geom_line(size=.25) +
  geom_point(size=.25) +
  geom_line() +
  geom_point() +
  facet_wrap(~tau, labeller = label_bquote(cols = tau == .(tau)))
gg
