# Fit gamlss models with skewed-t residuals
# launched with 
# R CMD BATCH --no-restore app/gamlss-fit.R app/log_fit.txt


# packages ----------------------------------------------------------------
library(tidyverse)
library(gamlss)
library(parallel)
mc_cores <- 3


# load data ---------------------------------------------------------------
data <- readRDS("app/data/data_pp.rds")
x_cols <- readRDS("app/data/x_cols.rds")

# stations we are concerned with
stns_id <- unique(data$station_id)

# placeholder for residuals and fitted moments
data <- data %>% mutate(res = as.numeric(NA),
                        res_mod = as.numeric(NA),
                        mu = as.numeric(NA),
                        sig = as.numeric(NA),
                        nu = as.numeric(NA),
                        tau = as.numeric(NA))


# fit gamlss model to each station ----------------------------------------

# formula (lots of interactions, but we have a lot of data)
ff <- temp ~ temp1*temp2*temp3 * (s1+c1+s2+c2+s3+c3+s4+c4)

dir.create("app/params") # folder for intermediate parameters

# Loop strategy
# for(stn in stns_id){
#   
#   cat("\n")
#   cat("--------------------------\n")
#   cat("--------------------------\n")
#   cat("*** considering station (id) ", stn, "\n \n")
#   
#   ind <- which(data$keep & data$station_id == stn)
# 
#   cat("fit model\n \n")
#   mod <- gamlss(formula = ff, ff, ff, ff,
#                 family = gamlss.dist::ST1(),
#                 data = data[ind, x_cols],
#                 control = gamlss.control(c.crit = .01, n.cyc = 500),
#                 i.control = glim.control(cc = .01))
#   
#   mu <- mod$mu.coefficients
#   sig <- mod$sigma.coefficients
#   nu <- mod$nu.coefficients
#   tau <- mod$tau.coefficients
#   
#   saveRDS(cbind(mu, sig, nu, tau), paste0("app/params/P_", stn, ".rds"))
#   
#   cat("compute of residuals\n \n")
#   data[ind,]$res_mod <- mod$residuals
#   data[ind,]$mu <- c(mod$mu.x %*% mu)
#   data[ind,]$sig <- exp(c(mod$sigma.x %*% sig))
#   data[ind,]$nu <- c(mod$nu.x %*% nu)
#   data[ind,]$tau <- exp(c(mod$tau.x %*% tau))
#   
#   data[ind,]$res <- gamlss.dist::pST1(data[ind,]$temp,
#                                       data[ind,]$mu,
#                                       data[ind,]$sig,
#                                       data[ind,]$nu,
#                                       data[ind,]$tau)
#   
#   cat("saving backup data\n \n")
#   saveRDS(data, "app/data/data_res_backup.rds")
# }
# 
# saveRDS(data, "app/data/data_res.rds")
# 


# parallel strategy
data <- mclapply(stns_id, \(stn){
  
  # cat("\n")
  # cat("--------------------------\n")
  # cat("--------------------------\n")
  # cat("*** considering station (id) ", stn, "\n \n")
  # 
  data <- data %>% filter(station_id == stn)
  ind <- which(data$keep)
  
  # cat("fit model\n \n")
  mod <- gamlss(formula = ff, ff, ff, ff,
                family = gamlss.dist::ST1(),
                data = data[ind, x_cols],
                control = gamlss.control(c.crit = .01, n.cyc = 500),
                i.control = glim.control(cc = .01))
  
  mu <- mod$mu.coefficients
  sig <- mod$sigma.coefficients
  nu <- mod$nu.coefficients
  tau <- mod$tau.coefficients
  
  saveRDS(cbind(mu, sig, nu, tau), paste0("app/params/P_", stn, ".rds"))
  
  # cat("compute of residuals\n \n")
  data[ind,]$res_mod <- mod$residuals
  data[ind,]$mu <- c(mod$mu.x %*% mu)
  data[ind,]$sig <- exp(c(mod$sigma.x %*% sig))
  data[ind,]$nu <- c(mod$nu.x %*% nu)
  data[ind,]$tau <- exp(c(mod$tau.x %*% tau))
  
  data[ind,]$res <- gamlss.dist::pST1(data[ind,]$temp,
                                      data[ind,]$mu,
                                      data[ind,]$sig,
                                      data[ind,]$nu,
                                      data[ind,]$tau)
  
  # cat("saving backup data\n \n")
  saveRDS(data, paste0("app/data/data_res_", stn, ".rds"))
  
  return(data)
}, mc.cores = mc_cores) %>% do.call(what = "rbind")


saveRDS(data, "app/data/data_res.rds")



# # if needed
# data <- lapply(stns_id, \(stn){
#   readRDS(paste0("app/data/data_res_", stn, ".rds"))
# }) %>% do.call(what = "rbind")
# saveRDS(data, "app/data/data_res.rds")
