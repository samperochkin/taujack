# This is meant to be run after line 237 of "app/2-data-analysis.R"

# focus on p=15
L_max <- 14
labs <- paste0("k=", 1:L_max + 1)
names(labs) <- 1:L_max + 1

# for 99% confidence int.
a <- abs(qnorm(.01,0,1))

gg <- data.frame(city = rep(stns_name, each=L),
                 k = rep(1:L + 1, times=q), p = p,
                 tau = c(Ths),
                 sigma = sqrt(c(sapply(1:q, \(k) diag(Shs[[k]])/ns[k])))) %>%
  filter(k <= L_max + 1) %>%
  ggplot(aes(x = city, y = tau, col=city, fill=city)) +
  theme_bw() + cc +
  ylab(bquote(hat(tau)["k"])) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = c(.75,.6),
        text = element_text(size=9),
        axis.text = element_text(size=9),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12, angle = 0, vjust = 0.5),
        legend.background = element_rect(colour="black"),
        legend.title = element_text(size = 9),
        legend.text=element_text(size=9)) +
  geom_errorbar(aes(ymin = tau - a*sigma, ymax = tau + a*sigma)) +
  facet_grid(~k, labeller = as_labeller(labs))
  # geom_ribbon(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma), alpha=.3, colour=NA) +
  # geom_line(size=.25) + geom_point(size=.25)
gg

# for paper (cc <- NULL)
ggsave(filename = "figures/tau_example.pdf", device = "pdf",
       width = 6, height = 2.5, units = "in")

# gg + theme(legend.position = "none") + coord_cartesian(xlim = c(1,3), ylim = c(.3,.5))

