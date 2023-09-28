library(lunar)
library(dplyr)

stns_id <- readRDS("app/data/stns_id.rds")
cols <- colorRampPalette(c("red", "yellow", "green"))(101)

for(i in seq_along(stns_id)){
  stn <- stns_id[i]
  data <- readRDS(paste0("app/data/daily_data_", stn, ".rds"))
  
  # time variables setup
  data <- data %>% mutate(date = as.Date(date),
                          year = as.integer(year),
                          month = as.integer(month),
                          day = as.integer(day))
  
  data <- data %>% mutate(time = as.integer(date),
                          yday = lubridate::yday(date),
                          week = lubridate::week(date),
                          seas = lunar::terrestrial.season(date))
  
  data <- data %>% mutate(temp = max_temp)
  
  data0 <- data %>% group_by(year, month) %>% summarise(prop = mean(!is.na(temp)))
  
  # data0 <- data0[1:200,]
  # plot(data0$year + data0$month/12, y = rep(i, nrow(data0)), col=cols[round(100*data0$prop)], cex=.1, pch=19)
  
  if(i == 1){
    plot(data0$year + data0$month/12, y = rep(i, nrow(data0)), col=cols[1+round(100*data0$prop)], cex=.1,
         xlim = c(1830, 2024), ylim = c(1, length(stns_id)))
  }else{
    points(data0$year + data0$month/12, y = rep(i, nrow(data0)), col=cols[1+round(100*data0$prop)], cex=.1)
  }
  text(x = 1830, y = i, labels = paste0(data$station_name[1], " (", data$station_id[1], ")"),
       pos = 4, cex = .6)
}


