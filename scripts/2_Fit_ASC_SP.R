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
    summarise(n_sp = sum(n_sp), Pop = sum(Pop)) %>%
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
    summarise(DeaR = weighted.mean(DeaR, Noti), pr_sp = 1) %>%
    arrange(Sex, Year)
  
  yrs <- sort(unique(noti$Year))
  n_t <- length(yrs)
  
  ## Asymptomatic -> Symptomatic -> Notification ----
  dat <- list(
    YearSurveyed = prv$Year[1],
    N = round(prv$N),
    Asym = round(prv$AsymSp),
    Sym = round(prv$SymSp),
    Years = yrs,
    Pop = t(matrix(noti$Pop, n_t)),
    Noti = t(matrix(noti$n_sp, n_t)),
    n_t = n_t,
    n_gp = nrow(prv)
  )
  
  exo <- get_exo_asc(mor, untr_a = "no_untr", untr_s = "avg", bg_death = T, pr_fp = 0)
  
  dat_sp <- c(dat, exo)
  
  model <- readRDS(file = "stan/m3_as_duration_free.rds")
  fitted_sp <- sampling(model, data = dat_sp, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  summary(fitted_sp, pars = c("r_sym", "r_det", "adr", "dur_a", "dur_s"))$summary
  
  save(fitted_sp, dat_sp, file = paste0("out/ASC/Post_Sp_", iso, ".rdata"))
}

