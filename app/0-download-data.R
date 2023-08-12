
###########################################################################
# Script dowloading data from Weather Canada ------------------------------
###########################################################################



# packages ----------------------------------------------------------------

# install.packages("weathercan",  repos = c("https://ropensci.r-universe.dev", "https://cloud.r-project.org"))
# install.packages(c('lutz', 'sf')) # dependencies
library(dplyr)
library(weathercan)



# select stations ---------------------------------------------------------

?stations
names(stations())
stations() %>%
  filter(start <= 1900, end >= 2000, interval == "day") %>%
  print(n=nrow(.))

# Interesting selections:

# airports: !is.na(TC_id)
# stns <- stations() %>%
#   filter(start <= 1900, end >= 2000, interval == "day", !is.na(TC_id))
# stns_id <- stns %>%
#   dplyr::select(station_id) %>%
#   unlist(use.names = F)

# Ottawa and Toronto:
# stns_id <- c(271, 568, 707, 2925, 4333, 4712, 5051)
stns_id <- c(271, 568, 707, 4333, 5051)
stns <- stations() %>%
  filter(station_id %in% stns_id, interval == "day")
stns


# download raw data -------------------------------------------------------

# from terminal:
# R CMD BATCH --no-restore app/download-data.R app/log_dd.txt

for(stn in stns_id){
  data <- weather_dl(stn, interval = "day", trim = F, verbose = T)
  saveRDS(data, paste0("app/data/daily_data_", stn, ".rds"))
}

# record successful downloads
stns_id <- sapply(strsplit(list.files("app/data/", "daily_data"), "_"), \(x) x[3])
stns_id <- as.numeric(sapply(strsplit(stns_id, "\\."), \(x) x[1]))
saveRDS(stns_id, "app/data/stns_id.rds")
