
# seas part originally appended to the data analysis.
# Decided to do a multivariate thing instead.
# Needs to be corrected if ever used


# Seasonal analysis -------------------------------------------------------
seass <- levels(data$seas)
d <- 30
K <- 25

res <- lapply(seass, \(s){
  X0 <- data %>% ungroup %>%
    select(date, station_id, pseudo_temp) %>%
    mutate(pseudo_temp_seas = ifelse(data$seas == s, pseudo_temp, NA)) %>%
    pivot_wider(names_from = station_id, values_from = pseudo_temp_seas) %>%
    select(any_of(as.character(stns_id))) %>%
    as.matrix()
  
  Y0 <- sapply(1:d, \(k) X0[k:(nrow(X)-d+k)])
  na_rows <- which(apply(is.na(Y0), 1, any))
  n0 <- nrow(Y0) - length(na_rows)
  
  C0 <- matrix(NA, nrow(Y0), d-1)
  C0[-na_rows,] <- dac_serial(Y0[-na_rows,])
  ths0 <- (2^(2:d-1) * colSums(C0, na.rm=T)/(n0*(n0-1)) - 1)/(2^(2:d-1) - 1)
  g0 <- t((2^(2:d-1) * t(C0)/(n0-1) - 1)/(2^(2:d-1) - 1) - ths0)
  
  # compute sigs
  zetas0 <- sapply(0:K, \(k){
    colMeans(g0[1:(n0-k),]*g0[(k+1):n0,], na.rm=T)
  }) |> t()
  plot(0:K, zetas0[,1], type="o", xlim=c(0,15), main = s)
  for(i in 2:(d-1)) lines(0:K, zetas0[,i], type="o", col=i)
  abline(h=0, lty=2)
  sigs0 <- zetas0[1,] + 2*colSums(zetas0[-1,])
  
  # sigma already normalized
  list(ths = ths0, sigs = sqrt(c(sigs0)/n0))
})

# check taus
thss <- sapply(res, "[[", "ths")
plot(thss[,1], ylim=c(0,max(thss)), type="o")
for(i in 2:4) lines(thss[,i], col=i, type="o")
abline(h=0, lty=2)

# check sigs
sigss <- sapply(res, "[[", "sigs")
plot(sigss[,1], ylim=c(0,max(sigss)), type="o")
for(i in 2:4) lines(sigss[,i], col=i, type="o")
abline(h=0, lty=2)

# plots
df <- data.frame(tau = c(thss),
                 sigma = c(sigss),
                 dim = rep(2:d, times=4),
                 seas = rep(seass, each=d-1))

ggplot(df, aes(x = dim, y = tau, col=seas, group=seas)) +
  theme_light() +
  geom_line() +
  # geom_errorbar(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma)) +
  # geom_line(aes(y = tau - 3*sigma), lty=2, alpha=.5) +
  # geom_line(aes(y = tau + 3*sigma), lty=2, alpha=.5) +
  geom_point()

ggplot(df, aes(x = dim, y = tau, col=seas, group=seas)) +
  theme_light() +
  # coord_cartesian(xlim = c(2,5), ylim = c(.2, .45)) +
  # coord_cartesian(xlim = c(6,10), ylim = c(.05, .2)) +
  # coord_cartesian(xlim = c(11,20), ylim = c(0, .05)) +
  coord_cartesian(xlim = c(21,30), ylim = c(0, .01)) +
  geom_line() +
  # geom_errorbar(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma)) +
  geom_line(aes(y = tau - 3*sigma), lty=2, alpha=.5) +
  geom_line(aes(y = tau + 3*sigma), lty=2, alpha=.5) +
  geom_point()

ggplot(df, aes(x = dim, y = tau, col=seas, group=seas)) +
  theme_light() +
  coord_cartesian(xlim = c(d-10,d), ylim = c(-.0001,.005)) +
  geom_line() +
  # geom_errorbar(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma)) +
  geom_line(aes(y = tau - 3*sigma), lty=2, alpha=.5) +
  geom_line(aes(y = tau + 3*sigma), lty=2, alpha=.5) +
  geom_point()

ggplot(df %>% filter(dim %in% seq(2,30,7)), aes(x = tau, y = seas, col=seas, group=seas)) +
  theme_light() +
  geom_errorbar(aes(xmin = tau - 3*sigma, xmax = tau + 3*sigma)) +
  # geom_line(aes(y = tau - 3*sigma), lty=2, alpha=.5) +
  # geom_line(aes(y = tau + 3*sigma), lty=2, alpha=.5) +
  geom_point() +
  facet_wrap(~dim, nrow=1, scales = "free_x")



d <- 5
K <- 15
months <- 1:12
res <- lapply(months, \(m){
  cat("month:", m, "\n")
  X0 <- data %>% ungroup %>%
    select(date, station_id, pseudo_temp) %>%
    mutate(pseudo_temp_seas = ifelse(data$month == m, pseudo_temp, NA)) %>%
    pivot_wider(names_from = station_id, values_from = pseudo_temp_seas) %>%
    select(any_of(as.character(stns_id))) %>%
    as.matrix()
  
  Y0 <- sapply(1:d, \(k) X0[k:(nrow(X)-d+k)])
  na_rows <- which(apply(is.na(Y0), 1, any))
  n00 <- nrow(Y0)
  n0 <- n00 - length(na_rows)
  
  C0 <- matrix(NA, nrow(Y0), d-1)
  C0[-na_rows,] <- dac_serial(Y0[-na_rows,])
  ths0 <- (2^(2:d-1) * colSums(C0, na.rm=T)/(n0*(n0-1)) - 1)/(2^(2:d-1) - 1)
  g0 <- t((2^(2:d-1) * t(C0)/(n0-1) - 1)/(2^(2:d-1) - 1) - ths0)
  
  # compute sigs
  zetas0 <- sapply(0:K, \(k){
    colMeans(g0[1:(n00-k),]*g0[(k+1):n00,], na.rm=T)
  }) |> t()
  plot(0:K, zetas0[,1], type="o", xlim=c(0,15), main = paste0("month ",m))
  for(i in 2:(d-1)) lines(0:K, zetas0[,i], type="o", col=i)
  abline(h=0, lty=2)
  sigs0 <- zetas0[1,] + 2*colSums(zetas0[-1,])
  
  # sigss scaled
  list(ths = ths0, sigs = sqrt(c(sigs0)/n0))
})

# check taus
thss <- sapply(res, "[[", "ths")
plot(thss[,1], ylim=c(0,max(thss)), type="o")
for(i in 2:12) lines(thss[,i], col=i, type="o")
abline(h=0, lty=2)

# check sigs
sigss <- sapply(res, "[[", "sigs")
plot(sigss[,1], ylim=c(0,max(sigss)), type="o")
for(i in 2:12) lines(sigss[,i], col=i, type="o")
abline(h=0, lty=2)

# plots
df <- data.frame(tau = c(thss),
                 sigma = c(sigss),
                 dim = rep(2:d, times=4),
                 month = rep(months, each=d-1))

ggplot(df, aes(x = dim, y = tau, col=as.factor(month), group=month)) +
  theme_light() +
  geom_line() +
  # geom_errorbar(aes(ymin = tau - 2*sigma, ymax = tau + 2*sigma)) +
  geom_line(aes(y = tau - 3*sigma), lty=2, alpha=.5) +
  geom_line(aes(y = tau + 3*sigma), lty=2, alpha=.5) +
  geom_point()

ggplot(df, aes(x = dim, y = tau, col=as.factor(month), group=month)) +
  theme_light() +
  geom_line() +
  geom_line(aes(y = tau - 3*sigma), lty=2, alpha=.5) +
  geom_line(aes(y = tau + 3*sigma), lty=2, alpha=.5) +
  geom_point()

ggplot(df, aes(x = tau, y = month, col=as.factor(month), group=month)) +
  theme_light() +
  geom_errorbar(aes(xmin = tau - 3*sigma, xmax = tau + 3*sigma)) +
  # geom_line(aes(y = tau - 3*sigma), lty=2, alpha=.5) +
  # geom_line(aes(y = tau + 3*sigma), lty=2, alpha=.5) +
  geom_point() +
  facet_wrap(~dim, nrow=1)
