
# packages ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# load data ---------------------------------------------------------------
data <- readRDS("app/data/data_pp.rds")


# Global analysis ---------------------------------------------------------

names(data)
data0 <- data %>% select(date, station_id, temp) %>%
  pivot_wider(names_from = station_id, values_from = temp) %>% 
  drop_na()

X <- data0 %>% ungroup %>%
  # X <- data00 %>% ungroup %>% filter(seas == "Winter") %>%
  select(any_of(as.character(stns_id))) %>%
  as.matrix()
