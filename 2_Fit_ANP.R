library(tidyverse)
library(rstan)

options(mc.cores = min(5, parallel::detectCores()))


source("R/var_exo_anp.R")
source("data/country_list.R")


###
n_iter = 5000
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
    summarise(n_sp = round(sum(n_sp)), 
              n_sn = round(sum(n_sn)), 
              Pop = round(sum(Pop))) %>%
    arrange(Sex)

  
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
    N = prv$N,
    Asym = prv$Asym,
    Sn = prv$SymSn,
    Sp = prv$SymSp,
    Years = yrs,
    Pop = t(matrix(noti$Pop, n_t)),
    NotiSn = t(matrix(noti$n_sn, n_t)),
    NotiSp = t(matrix(noti$n_sp, n_t)),
    n_t = n_t,
    n_gp = nrow(prv)
  )
  
  exo <- get_exo_anp(mor, untr_a = "no_untr", untr_s = "full", bg_death = T, pr_fp = 0)
  
  dat_anp <- c(dat, exo)
  
  
  model <- readRDS(file = "stan/m1_duration_free.rds")
  fitted_anp_uni <- sampling(model, data = dat_anp, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  check_divergences(fitted_anp_uni)
  summary(fitted_anp_uni, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary
  
  pars <- as_tibble(extract(fitted_anp_uni, c("p_sp", "r_tr")))  
  plot(pars$r_tr, pars$p_sp, main = iso)
  
  model <- readRDS(file = "stan/m1_reg.rds")
  fitted_anp_reg <- sampling(model, data = dat_anp, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  check_divergences(fitted_anp_reg)
  summary(fitted_anp_reg, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr", "lrr_cs_sn", "lrr_cs_sp", "lrr_sym"))$summary
  
  
  save(fitted_anp_uni, fitted_anp_reg, dat_anp, file = paste0("out/ANP/Post_", iso, ".rdata"))
}



