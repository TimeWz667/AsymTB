library(tidyverse)
library(rstan)

options(mc.cores = min(5, parallel::detectCores()))


source("R/var_exo_anp.R")

###
n_iter = 5000
n_collect = 1000
n_chain = 3


model <- readRDS(file = "stan/m1_duration_free.rds")


### Blantyre
load("data/Input_BLT.rdata")


notification <- notification %>% filter(Year > 2014)

prv <- prevalence %>% 
  group_by(Year, HIV) %>%
  summarise(N = round(sum(N)), Asym = round(sum(Asym)), SymSn = round(sum(SymSn)), SymSp = round(sum(SymSp))) %>%
  arrange(HIV)


noti <- notification %>%
  group_by(Year, HIV) %>%
  summarise(n_all = round(sum(n_all)), 
            n_sp = round(sum(n_sp)), 
            n_sn = round(sum(n_sn)), 
            Pop = round(sum(Pop))) %>%
  arrange(HIV)


yrs <- sort(unique(noti$Year))
n_t <- length(yrs)

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


mor <- mortality %>% 
  left_join(tibble(Year = 2019, HIV = c("HIV", "NonHIV"))) %>%
  mutate(DeaR = ifelse(HIV == "HIV", DeaR + 11 / 1000, DeaR)) %>%
  group_by(Year, HIV) %>%
  summarise(DeaR = weighted.mean(DeaR, Noti)) %>%
  arrange(Year, HIV)


exos <- list(
  H1 = get_exo_anp_hiv(mor, untr_a = "no_untr", untr_s = "no_untr", untr_hiv = "as_nonhiv", bg_death = T, pr_fp = 0),
  H2 = get_exo_anp_hiv(mor, untr_a = "no_untr", untr_s = "full", untr_hiv = "as_nonhiv", bg_death = T, pr_fp = 0),
  H3 = get_exo_anp_hiv(mor, untr_a = "no_untr", untr_s = "full", untr_hiv = "full", bg_death = T, pr_fp = 0)
)


fitted_blt <- lapply(exos, function(exo) {
  fitted_hiv <- sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  print(summary(fitted_hiv, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary)
  fitted_hiv
})



#### Model fitting -----
### Kenya
load("data/Input_KEN_Full.rdata")


prv <- prevalence %>% 
  group_by(Year, HIV) %>%
  summarise(N = round(sum(N)), Asym = round(sum(Asym)), SymSn = round(sum(SymSn)), SymSp = round(sum(SymSp))) %>%
  arrange(HIV)


noti <- notification %>%
  group_by(Year, HIV) %>%
  summarise(n_all = round(sum(n_all)), 
            n_sp = round(sum(n_sp)), 
            n_sn = round(sum(n_sn)), 
            Pop = round(sum(Pop))) %>%
  arrange(HIV)


yrs <- sort(unique(noti$Year))
n_t <- length(yrs)

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

mor <- mortality %>% 
  left_join(tibble(Year = 2016, HIV = c("HIV", "NonHIV"))) %>%
  mutate(DeaR = ifelse(HIV == "HIV", DeaR + 21/1300, DeaR)) %>%
  group_by(Year, HIV) %>%
  summarise(DeaR = weighted.mean(DeaR, Noti)) %>%
  arrange(Year, HIV)


exos <- list(
  H1 = get_exo_anp_hiv(mor, untr_a = "no_untr", untr_s = "no_untr", untr_hiv = "as_nonhiv", bg_death = T, pr_fp = 0),
  H2 = get_exo_anp_hiv(mor, untr_a = "no_untr", untr_s = "full", untr_hiv = "as_nonhiv", bg_death = T, pr_fp = 0),
  H3 = get_exo_anp_hiv(mor, untr_a = "no_untr", untr_s = "full", untr_hiv = "full", bg_death = T, pr_fp = 0)
)


fitted_ken <- lapply(exos, function(exo) {
  fitted_hiv <- sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  print(summary(fitted_hiv, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary)
  fitted_hiv
})


save(fitted_blt, fitted_ken, file = "out/Sens/Post_HIV.rdata")


