library(tidyverse)
library(rstan)

options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = "-march=corei7 -mtune=corei7")

source("R/get_exo.R")


countries <- c(
  KHM = "Cambodia",
  KEN = "Kenya",
  LAO = "Lao People's Democratic Republic", 
  MWI = "Malawi", 
  PAK = "Pakistan", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda", #
  VNM = "Viet Nam", 
  ZMB = "Zambia"
)

countries_cs <- c(
  "Malawi", "Philippines", 
  "United Republic of Tanzania", "Uganda" , "Zambia"
)



for (i in 1:length(countries)) {
  iso <- names(countries)[i]
  country <- countries[i]
  
  print(country)
  
  load(paste0("data/Input_", iso, ".rdata"))
  
  prv <- prevalence
  noti <- notification %>% arrange(Sex, Year)
  dr <- mortality$DeaR
  
  yrs <- sort(unique(noti$Year))
  n_t <- length(yrs)

  
  
  ## Asymptomatic -> Symptomatic -> Notification ----
  dat <- list(
    YearSurveyed = prv$Year[1],
    N = round(prv$N),
    Asym = round(prv$Asym),
    Sym = round(prv$Sym),
    Years = yrs,
    Pop = t(matrix(noti$Pop, n_t)),
    Noti = t(matrix(noti$n_all, n_t)),
    n_t = n_t,
    n_gp = nrow(prv)
  )
  
  
  exo <- get_exo_sym(dr, "avg")
  
  dat_as <- c(dat, exo)
  
  model <- readRDS(file = "stan/m3_as_uni.rds")
  fitted_as_uni <- sampling(model, data = dat_as, iter = 3000, chain = 2)
  summary(fitted_as_uni, pars = c("r_sym", "r_det", "adr", "dur_a", "del_s"))$summary
  
  model <- readRDS(file = "stan/m3_as_fixed.rds")
  fitted_as_fixed <- sampling(model, data = dat_as, iter = 3000, chain = 2)
  summary(fitted_as_fixed, pars = c("r_sym", "r_det", "adr", "dur_a", "del_s"))$summary
  
  if (country %in% countries_cs) {
    ## Asymptomatic -> Symptomatic -> Care-seeking attempt -> Notification
    dat <- list(
      YearSurveyed = prv$Year[1],
      N = round(prv$N),
      Asym = round(prv$Asym),
      Sym = round(prv$SymNC),
      CS = round(prv$SymCS),
      Years = yrs,
      Pop = t(matrix(noti$Pop, n_t)),
      Noti = t(matrix(noti$n_all, n_t)),
      n_t = n_t,
      n_gp = nrow(prv)
    )
    
    exo <- get_exo_sym(dr, "avg")
    
    dat_asc <- c(dat, exo)
    
    model <- readRDS(file = "stan/m3_asc_uni.rds")
    fitted_asc_uni <- sampling(model, data = dat_asc, iter = 3000, chain = 2)
    summary(fitted_asc_uni, pars = c("r_sym", "r_aware", "r_det", "adr", "dur_a", "dur_s", "del_s", "del_c", "del"))$summary
    
    model <- readRDS(file = "stan/m3_asc_fixed.rds")
    fitted_asc_fixed <- sampling(model, data = dat_asc, iter = 3000, chain = 2)
    summary(fitted_asc_fixed, pars = c("r_sym", "r_aware", "r_det", "adr", "dur_a", "dur_s", "del_s", "del_c", "del"))$summary
    
    save(fitted_as_uni, fitted_as_fixed, dat_as,
         fitted_asc_uni, fitted_asc_fixed, dat_asc,
         file = paste0("output/ASC/Post_", iso, ".rdata"))
  } else {
    save(fitted_as_uni, fitted_as_fixed, dat_as,
         file = paste0("output/ASC/Post_", iso, ".rdata"))
  }
  
}



