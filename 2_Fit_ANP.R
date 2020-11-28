library(tidyverse)
library(rstan)

options(mc.cores = min(5, parallel::detectCores()))


source("R/get_exo.R")


###
n_iter = 5000
n_collect = 1000
n_chain = 3

###
source("data/country_list.R")



for (i in 1:length(countries)) {
  iso <- glue::as_glue(names(countries)[i])
  country <- countries[i]
  
  print(country)
  
  load("data/Input_" + iso + ".rdata")
  
  prv <- prevalence %>% 
    group_by(Year, Sex) %>%
    summarise(N = round(sum(N)), Asym = round(sum(Asym)), SymSn = round(sum(SymSn)), SymSp = round(sum(SymSp))) %>%
    arrange(Sex)
  
  
  noti <- notification %>%
    group_by(Year, Sex) %>%
    summarise(n_all = round(sum(n_all)), 
              n_sp = round(sum(n_sp)), 
              n_sn = round(sum(n_sn)), 
              Pop = round(sum(Pop))) %>%
    arrange(Sex)

  
  if (iso == "KHM") {
    noti <- noti %>% filter(Year >= 2014)
  }
  if (iso == "UGA") {
    noti <- noti %>% filter(Year >= 2014)
  }
  if (iso == "VNM") {
    noti <- noti %>% filter(Year >= 2016)
  }
  
  
  dr <- mortality %>%
    group_by(Year, Sex) %>%
    summarise(DeaR = sum(DeaR * Pop) / sum(Pop)) %>%
    arrange(Sex)
  dr <- dr$DeaR

  
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
  
  exo <- get_exo(rep(F, nrow(prv)), dr, untr = "full", bg_death = T)
  
  dat_anp <- c(dat, exo)
  
  
  model <- readRDS(file = "stan/m1_duration_free.rds")
  fitted_anp_uni <- sampling(model, data = dat_anp, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  check_divergences(fitted_anp_uni)
  summary(fitted_anp_uni, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary
  
  pars <- as_tibble(extract(fitted_anp_uni, c("p_sp", "r_tr")))  
  plot(pars$p_sp, pars$r_tr, main = iso)
  
  model <- readRDS(file = "stan/m1_reg.rds")
  fitted_anp_reg <- sampling(model, data = dat_anp, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  check_divergences(fitted_anp_reg)
  summary(fitted_anp_reg, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr", "lrr_cs_sn", "lrr_cs_sp", "lrr_sym"))$summary
  
  
  save(fitted_anp_uni, fitted_anp_reg, dat_anp, file = paste0("out/ANP/Post_", iso, ".rdata"))
}



