
###########################################################################
# Script preprocessing the data -------------------------------------------
###########################################################################

# launched with 
# R CMD BATCH --vanilla --no-restore app/1-preprocessing.R app/log_pp.txt

# packages ----------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(ggplot2)
library(spatstat.geom) # for weighted ecdf




# load data and minor changes to names ------------------------------------

# stations we consider
stns_id <- readRDS("app/data/stns_id.rds")
p <- length(stns_id)

# meta info about the stations
meta <- weathercan::stations() %>%
  filter(station_id %in% stns_id, interval == "day") %>%
  print(nrow = nrow(.))
stns_name <- sapply(stns_id, \(i) meta %>% filter(station_id == i) %>% select(station_name) %>% unlist())

# load data
filepaths <- paste0("app/data/daily_data_", stns_id, ".rds")
data <- lapply(filepaths, readRDS) |> dplyr::bind_rows()

# clean station names
data[data$station_name == "OTTAWA CDA",]$station_name <- "Ottawa"
data[data$station_name == "WELLAND",]$station_name <- "Welland"
data[data$station_name == "TORONTO",]$station_name <- "Toronto"
stns_name[stns_name == "OTTAWA CDA"] <- "Ottawa"
stns_name[stns_name == "WELLAND"] <- "Welland"
stns_name[stns_name == "TORONTO"] <- "Toronto"
str(data)



# setup time variables ----------------------------------------------------
data <- data %>% 
  select(station_name, station_id, date, mean_temp) %>%
  mutate(date = as.Date(date))

# use same timespan for all stations
dates <- seq(min(data$date),max(data$date), by=1)
dates <- tibble(date = rep(dates, length(stns_id)),
                station_id = rep(stns_id, each = length(dates)))
data <- right_join(data, dates, by=c("date", "station_id"))
rm(dates)

# variable creation
data <- data %>% 
  mutate(time = as.integer(date), year = lubridate::year(date),
         month = lubridate::month(date), week = lubridate::week(date),
         day = lubridate::mday(date), yday = lubridate::yday(date),
         station_name = stns_name[match(station_id, stns_id)]) %>%
  arrange(station_id, date)



# summary plots and final choice of dates considered ----------------------

# clear time trends (w/ NAs)
data %>% group_by(station_name, year) %>%
  summarise(mm_temp = mean(mean_temp, na.rm=T)) %>%
  ggplot(aes(x=year, y=mm_temp)) +
  theme_bw() + ylab("(yearly) mean temperature") +
  geom_line(size=.25) + geom_point(size=.25) +
  facet_grid(rows="station_name", scales = "free_y")
# clear time trends (w/o NAs)
gg <- data %>% group_by(station_name, year) %>%
  summarise(mm_temp = mean(mean_temp, na.rm=F)) %>%
  ggplot(aes(x=year, y=mm_temp)) +
  theme_bw() + ylab("(yearly) mean temperature") +
  geom_line(size=.25) + geom_point(size=.25) +
  facet_grid(rows="station_name", scales = "free_y")
gg # will reuse gg later
ggsave(filename = "app/figures/series_raw_mean.pdf", device = "pdf",
       width = 6.5, height = 4, units = "in")

# clear seasonal trends
data %>% filter(year %in% 1975:1980) %>%
  group_by(station_name, month, year) %>%
  summarise(mm_temp = mean(mean_temp, na.rm=F)) %>%
  ggplot(aes(x=year+month/12, y=mm_temp)) +
  theme_bw() + ylab("(monthly) mean temperature") +
  geom_line(size=.25) + geom_point(size=.25) +
  facet_grid(rows="station_name", scales = "free_y")
ggsave(filename = "app/figures/series_raw_mean_seas.pdf", device = "pdf",
       width = 6.5, height = 4, units = "in")




# model mean based on time (seas. and long term trends) -------------------
year_len <- 365 + 6/24 + 9/60/24 + 9/60^2/24

# create sines-cosines basis (for seasonal trend)
nSC <- 6
s_cols <- paste0("s", 1:nSC)
c_cols <- paste0("c", 1:nSC)
for (i in 1:nSC) data <- data %>%
  mutate(!!s_cols[i] := sin(i*2*pi*time/year_len),
         !!c_cols[i] := cos(i*2*pi*time/year_len))

# create natural spline basis (for long term time trend) based on gg
gg
knot_years <- c(1882, 1925, 1953, 1997)
gg + geom_vline(xintercept = knot_years, lty=2)

degree <- 3
nB <- length(knot_years) + degree
b_cols <- paste0("b", 1:nB)
for(b in b_cols) data <- data %>% mutate(!!b := NA)

b_colss <- list()
for(k in 1:p){
  data0 <- data %>% filter(station_id == stns_id[k])
  min_time <- data0 %>% filter(!is.na(mean_temp)) %>% select(time) %>% unlist() %>% min()
  knot_times <- data0 %>% filter(day == 1, month == 1, year %in% knot_years, time >= min_time) %>%
    select(time) %>% unlist()
  l <- length(knot_times) + degree
  b_colss[[k]] <- b_cols[1:l]
  data[data$station_id == stns_id[k], b_colss[[k]]] <- splines::bs(data0$time, knots = knot_times, degree = degree)
}
gc()

# fit least-squares
# non_na <- !is.na(data$mean_temp) # for later
sc_vars <- paste0("(", paste0(c(s_cols, c_cols), collapse = "+"), ")")
lms <- lapply(1:p, \(k){
  
  all_cols <- c("station_id", "mean_temp", s_cols, c_cols, b_colss[[k]]) 
  b_vars <- paste0("(", paste0(b_colss[[k]], collapse = "+"), ")")
  
  print(k)
  # if(stns_name[k] == "Welland"){
  #   b_vars2 <- paste0("(", paste0("b", 2:nB, collapse = "+"),")")
  #   ff <- paste0("mean_temp ~ ", sc_vars, " + ", b_vars2, " + (s1+c1):", b_vars2)
  # }else if(stns_name[k] == "Toronto"){
  #   ff <- paste0("mean_temp ~ ", sc_vars, " + ", b_vars, " + (s1+c1):(", paste0("b", 1:(nB-1), collapse = "+"), ")")
  # }else if(stns_name[k] == "Ottawa"){
  #   b_vars2 <- paste0("(", paste0("b", 1:(nB-1), collapse = "+"),")")
  #   ff <- paste0("mean_temp ~ ", sc_vars, " + ", b_vars, " + (s1+c1):", b_vars)
  # }
  ff <- paste0("mean_temp ~ ", sc_vars, " + ", b_vars, " + (s1+c1):", b_vars)
  lm(formula = as.formula(ff),
     data = data[,all_cols] %>% filter(station_id == stns_id[k], !is.na(mean_temp)))
})
lapply(lms, \(l) summary(l))

# register the results
data <- data %>% mutate(mu = NA , ctemp = NA)
for(k in 1:p){
  ind <- data$station_id == stns_id[k] & !is.na(data$mean_temp)
  data[ind,] <- data[ind,] %>%  mutate(mu = lms[[k]]$fitted.values, ctemp = mean_temp - mu)
}

# remove b cols of data (for memory)
data <- data %>% 
  select(!any_of(b_cols[[which.max(sapply(b_cols, length))]])) %>%
  select(!any_of(c(s_cols, c_cols)))
gc()

# quick check (yearly)
data %>% group_by(year, station_name) %>%
  summarise(m_mt = mean(mean_temp, na.rm=F), m_mu = mean(mu, na.rm=F)) %>%
  # summarise(m_mt = mean(mean_temp, na.rm=T), m_mu = mean(mu, na.rm=T)) %>%
  ggplot(aes(x=year, y=m_mt)) +
  theme_bw() + theme(legend.position = "none") + ylab("(yearly) mean temperature") +
  geom_line(aes(y=m_mu), alpha=.35, size=.5, col="red") +
  geom_line(size=.25) + geom_point(size=.25) +
  facet_grid(rows = "station_name", scales = "free_y")
ggsave(filename = "app/figures/series_fit_yearly.pdf", device = "pdf",
       width = 6.5, height = 3, units = "in")

# quick check (monthly)
y0 <- 1945
data %>% filter(year %in% (y0 + 1:10)) %>%
  group_by(year, month, station_name) %>%
  summarise(m_mt = mean(mean_temp, na.rm=F),
            m_mu = mean(mu, na.rm=F)) %>%
  ggplot(aes(x=year+month/12, y=m_mt)) +
  theme_bw() + theme(legend.position = "none") + ylab("(seasonal) mean temperature") +
  geom_line(aes(y=m_mu), alpha=.35, size=.5, col="red") + geom_point(size=.25) +
  facet_grid(rows = "station_name", scales = "free_y")
ggsave(filename = "app/figures/series_fit_seas.pdf", device = "pdf",
       width = 6.5, height = 3, units = "in")

# quick check (monthly)
y0 <- 1950:1953
data %>% filter(year %in% y0) %>%
  ggplot(aes(x=date, y=mean_temp)) +
  theme_bw() + theme(legend.position = "none") + 
  ylab("centered temperature") + xlab("day of the year") +
  geom_line(aes(y=mu), alpha=.75, size=.5, col="red") +
  geom_point(size=.1) +
  facet_grid(rows = "station_name", scales = "free_y")
ggsave(filename = "app/figures/temp_daily.pdf", device = "pdf",
       width = 6.5, height = 3, units = "in")

# quick check (by yday)
y0 <- 1945 + 0:2
data %>% filter(year %in% y0) %>%
  ggplot(aes(x=date, y=ctemp)) +
  theme_bw() + theme(legend.position = "none") + 
  ylab("centered temperature") + xlab("day of the year") +
  geom_point(size=.1) +
  facet_grid(rows = "station_name", scales = "free_y")
ggsave(filename = "app/figures/ctemp_daily.pdf", device = "pdf",
       width = 6.5, height = 3, units = "in")
ggplot(data, aes(x=yday, y=ctemp)) +
  theme_bw() + theme(legend.position = "none") + 
  ylab("centered temperature") + xlab("day of the year") +
  geom_point(size=.1) +
  facet_grid(rows = "station_name", scales = "free_y")

# subsample for figure in paper
set.seed(34)
ggplot(data[sample(nrow(data), nrow(data)/5),], aes(x=yday, y=ctemp)) +
  theme_bw() + theme(legend.position = "none") + 
  ylab("centered temperature") + xlab("day of the year") +
  geom_point(size=.1) +
  facet_grid(rows = "station_name", scales = "free_y")
ggsave(filename = "app/figures/ctemp_seas.pdf", device = "pdf",
       width = 6.5, height = 3, units = "in")



# Compute pseudo-observations using ECDF ----------------------------------
data <- data %>% mutate(pseudo_temp = NA, diff = NA, w1 = NA, w = NA)
for(k in 1:p){
  id <- stns_id[k]
  cat("Working on ECDF for station ", id, ".\n")
  row_id <- data$station_id == id & !is.na(data$ctemp)
  
  for(tt in 1:366){
    if(tt %% 20 == 0) cat("progess: ", round(tt/366*100,2), "%\n")
    
    # construct weights based on yday using gaussian kernel
    data[row_id,] <- data[row_id,] %>%
      mutate(diff = tt - (time %% year_len)) %>%
      mutate(diff = abs(pmin(diff, year_len - diff))) %>%
      mutate(w1 = dnorm(diff,0,2)) %>%
      mutate(w = w1/sum(w1, na.rm = T)) 
    
    # construct ecdf function using all the data
    wecdf <- spatstat.geom::ewcdf(data[row_id,]$ctemp, data[row_id,]$w)
    
    # compute ecdf only for yday we care for
    ind <- which(data[row_id,]$yday == tt & !is.na(data[row_id,]$ctemp))
    data[row_id,][ind,] <- data[row_id,][ind,] %>%
      mutate(pseudo_temp = wecdf(ctemp))
    
    # heuristic correction for the largest (pseudo_obs = 1)
    # to avoid modifying all values (say, with pseudo_obs = pseudo_obs*n_ind/(n_ind+1))
    n_ind <- length(ind)
    data[row_id,][ind,] <- data[row_id,][ind,] %>%
      mutate(pseudo_temp = ifelse(pseudo_temp == 1, 1-1/(2*n_ind), pseudo_temp))
  }
}

# resulting histograms are even better than expected
ggplot(data, aes(x = pseudo_temp, fill=as.factor(yday))) +
  theme_bw() + theme(legend.position = "none") + xlab("(preprocessed) temperature") +
  geom_histogram(breaks=seq(0,1,.05), position = position_stack()) +
  facet_wrap(~station_name)
ggplot(data, aes(x = pseudo_temp)) +
  theme_bw() + theme(legend.position = "none") + xlab("(preprocessed) temperature") +
  scale_x_continuous(breaks = seq(0,1,.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  geom_histogram(breaks=seq(0,1,.05), position = position_stack()) +
  facet_wrap(~station_name)
ggsave(filename = "app/figures/pseudo_hist.pdf", device = "pdf",
       width = 6.5, height = 2, units = "in")

# note that there are equalities, but very very few of them. can disregard.
sapply(1:p, \(k){
  r <- data$station_id == id & !is.na(data$ctemp)
  tab <- table(data[r,]$pseudo_temp)
  tab[tab > 1]
})



# Save for later use ------------------------------------------------------
saveRDS(data, "app/data/data_pp.rds")
