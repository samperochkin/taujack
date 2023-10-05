library(ggplot2)
library(tidyverse)

n <- 20
set.seed(15)
data <- data.frame(x=sample(n)/(n+1), y=(1:n)/(n+1)) %>%
  mutate(gamma = x <= .5)

scoreFun <- function(gamma){
  n <- length(gamma)
  sapply(1:n, \(k){
    s1 <- sum(gamma[1:k])
    s2 <- sum(c(gamma[-(1:k)],0))
    s1*(n-k-s2) + (k-s1)*s2
  })
}
s <- which.max(scoreFun(data$gamma))/(n+1) + .5/(n+1)

df <- data.frame(x = c(.05,.95), y = c(.925,.925), lab = c("I[0]", "I[1]"))
ggplot(df, aes(x=x,y=y)) + 
  theme_bw() +
  theme(axis.title = element_text(size=10),
        axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), aspect.ratio = 1) +
  xlab(bquote(X[1])) + ylab(bquote(X[2])) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  geom_point(data=data, aes(x=x, y=y), col="gray", size=.5) +
  geom_vline(xintercept=.5, lty=2, size=.25) +
  geom_text(aes(label=lab), parse = T, size = 3.25)
ggsave("figures/dac-toy-1.pdf", 
       device = "pdf", height = 1.75, width = 1.75)


df <- data.frame(
  x = c(.05,.95, .05,.95),
  y = c(.075,.075, .925, .925),
  lab = c("I[0]^0", "I[1]^0", "I[0]^1", "I[1]^1")
)

ggplot(df, aes(x=x,y=y)) + 
  theme_bw() +
  theme(axis.title = element_text(size=10),
        axis.text.y = element_text(size=8.5, colour="black"), axis.ticks.y = element_line(size = .35),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(), aspect.ratio = 1) +
  xlab(bquote(X["k-1"])) + ylab(bquote(X["k"])) +
  scale_x_continuous(limits = c(0,1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,1), breaks = s, labels = "s", expand = c(0, 0)) +
  geom_point(data=data, aes(x=x, y=y), col="gray", size=.5) +
  geom_hline(yintercept=s, lty=2, size=.25) +
  geom_vline(xintercept=.5, lty=2, size=.25) +
  geom_text(aes(label=lab), parse = T, size = 3.25)
ggsave("figures/dac-toy-2.pdf", 
       device = "pdf", height = 1.75, width = 1.825)




