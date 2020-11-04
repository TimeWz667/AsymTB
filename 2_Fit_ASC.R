library(tidyverse)
library(rstan)

options(mc.cores = min(5, parallel::detectCores()))


source("R/get_exo.R")


###
n_iter = 3000
n_collect = 1000
n_chain = 3

###


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
  KEN = "Kenya",
  MWI = "Malawi", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda", #
  ZMB = "Zambia"
)



for (i in 1:length(countries)) {
  iso <- glue::as_glue(names(countries)[i])
  country <- countries[i]
  
  print(country)
  
  load("data/Input_" + iso + ".rdata")
  
  if (country %in% countries_cs) {
    prv <- prevalence %>%
      group_by(Year, Sex) %>%
      summarise(N = sum(N),
                Asym = sum(Asym),
                Sym = sum(Sym),
                SymNC = sum(SymNC),
                SymCS = sum(SymCS))
  } else {
    prv <- prevalence %>%
      group_by(Year, Sex) %>%
      summarise(N = sum(N),
                Asym = sum(Asym),
                Sym = sum(Sym))
  }
  
  
  noti <- notification %>% 
    group_by(Year, Sex) %>%
    summarise(n_all = sum(n_all), Pop = sum(Pop)) %>%
    arrange(Sex, Year)
    
  
  if (iso == "KHM") {
    noti <- noti %>% filter(Year >= 2014)
  }
  if (iso == "UGA") {
    noti <- noti %>% filter(Year >= 2014)
  }
  if (iso == "VNM") {
    noti <- noti %>% filter(Year >= 2016)
  }
  
  dr <- (mortality %>%
           group_by(Year, Sex) %>%
           summarise(DeaR = sum(DeaR * Pop) / sum(Pop), Pop = sum(Pop)))$DeaR
  
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
  fitted_as_uni <- sampling(model, data = dat_as, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  summary(fitted_as_uni, pars = c("r_sym", "r_det", "adr", "dur_a", "dur_s"))$summary
  check_divergences(fitted_as_uni)
  
  
  model <- readRDS(file = "stan/m3_as_fixed.rds")
  fitted_as_fixed <- sampling(model, data = dat_as, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  check_divergences(fitted_as_fixed)
  summary(fitted_as_fixed, pars = c("r_sym", "r_det", "adr", "dur_a", "dur_s"))$summary
  
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
    fitted_asc_uni <- sampling(model, data = dat_asc, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
    summary(fitted_asc_uni, pars = c("r_sym", "r_aware", "r_det", "adr", "dur_a", "dur_s", "dur_c"))$summary
    
    model <- readRDS(file = "stan/m3_asc_fixed.rds")
    fitted_asc_fixed <- sampling(model, data = dat_asc, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
    summary(fitted_asc_fixed, pars = c("r_sym", "r_aware", "r_det", "adr", "dur_a", "dur_s", "dur_c"))$summary
    
    save(fitted_as_uni, fitted_as_fixed, dat_as,
         fitted_asc_uni, fitted_asc_fixed, dat_asc,
         file = paste0("out/ASC/Post_", iso, ".rdata"))
  } else {
    save(fitted_as_uni, fitted_as_fixed, dat_as,
         file = paste0("out/ASC/Post_", iso, ".rdata"))
  }
  
}



