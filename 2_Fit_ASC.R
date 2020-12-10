library(tidyverse)
library(rstan)

options(mc.cores = min(5, parallel::detectCores()))


source("R/var_exo_asc.R")
source("data/country_list.R")


###
n_iter = 3000
n_collect = 1000
n_chain = 3


###
for (i in 1:length(countries)) {
  iso <- glue::as_glue(names(countries)[i])
  country <- countries[i]
  
  print(country)
  
  load("data/Input_" + iso + ".rdata")
  
  prv <- prevalence %>% arrange(Sex)
  
  noti <- notification %>% 
    group_by(Year, Sex) %>%
    summarise(n_all = sum(n_all), Pop = sum(Pop)) %>%
    arrange(Sex, Year)
    
  
  if (iso == "KHM") {
    noti <- noti %>% filter(Year >= 2014)
  } else if (iso == "UGA") {
    noti <- noti %>% filter(Year >= 2014)
  } else if (iso == "VNM") {
    noti <- noti %>% filter(Year >= 2016)
  } else if (iso == "PHL") {
    noti <- noti %>% filter(Year >= 2015)
  }
  
  
  mor <- mortality %>%
    group_by(Year, Sex) %>%
    summarise(DeaR = weighted.mean(DeaR, Noti), pr_sp = mean(pr_sp)) %>%
    arrange(Sex, Year)
    
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
  
  exo <- get_exo_asc(mor, untr_a = "no_untr", untr_s = "avg", bg_death = T, pr_fp = 0)
  
  dat_as <- c(dat, exo)
  
  model <- readRDS(file = "stan/m3_as_duration_free.rds")
  fitted_as_uni <- sampling(model, data = dat_as, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  summary(fitted_as_uni, pars = c("r_sym", "r_det", "adr", "dur_a", "dur_s"))$summary
  check_divergences(fitted_as_uni)
  
  
  #model <- readRDS(file = "stan/m3_as_fixed.rds")
  #fitted_as_fixed <- sampling(model, data = dat_as, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  #check_divergences(fitted_as_fixed)
  #summary(fitted_as_fixed, pars = c("r_sym", "r_det", "adr", "dur_a", "dur_s"))$summary
  
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
    
    exo <- get_exo_asc(mor, untr_a = "no_untr", untr_s = "avg", bg_death = T, pr_fp = 0)
    
    dat_asc <- c(dat, exo)
    
    model <- readRDS(file = "stan/m3_asc_duration_free.rds")
    fitted_asc_uni <- sampling(model, data = dat_asc, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
    summary(fitted_asc_uni, pars = c("r_sym", "r_aware", "r_det", "adr", "dur_a", "dur_s", "dur_c"))$summary
    check_divergences(fitted_asc_uni)
    
    #model <- readRDS(file = "stan/m3_asc_fixed.rds")
    #fitted_asc_fixed <- sampling(model, data = dat_asc, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
    #summary(fitted_asc_fixed, pars = c("r_sym", "r_aware", "r_det", "adr", "dur_a", "dur_s", "dur_c"))$summary
    #check_divergences(fitted_asc_fixed)
    
    save(fitted_as_uni, dat_as,
         fitted_asc_uni, dat_asc,
         file = paste0("out/ASC/Post_", iso, ".rdata"))
  } else {
    save(fitted_as_uni, dat_as,
         file = paste0("out/ASC/Post_", iso, ".rdata"))
  }
}
