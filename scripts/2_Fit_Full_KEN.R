library(tidyverse)
library(rstan)

options(mc.cores = min(5, parallel::detectCores()))


source("R/var_exo_anp.R")
load("data/Input_KEN_Full.rdata")


###
n_iter = 5000
n_collect = 1000
n_chain = 3


#### Aggregated ----
prv <- prevalence

noti <- notification %>%
  group_by(Year) %>%
  summarise(n_all = round(sum(n_all)), 
            n_sp = round(sum(n_sp)), 
            n_sn = round(sum(n_sn)), 
            Pop = round(sum(Pop)))


yrs <- sort(unique(noti$Year))
n_t <- length(yrs)

dat <- list(
  YearSurveyed = prv$Year[1],
  N = round(sum(prv$N)),
  Asym = round(sum(prv$Asym)),
  Sn = round(sum(prv$SymSn)),
  Sp = round(sum(prv$SymSp)),
  Years = yrs,
  Pop = noti$Pop,
  NotiSn = noti$n_sn,
  NotiSp = noti$n_sp,
  n_t = n_t
)


mor <- mortality %>%
  group_by(Year) %>%
  summarise(DeaR = weighted.mean(DeaR, Noti)) %>%
  arrange(Year)


exo <- get_exo_anp(mor, untr_a = "no_untr", untr_s = "full", bg_death = T, pr_fp = 0)



model <- readRDS(file = "stan/m0.rds")
fitted_total <- sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
check_divergences(fitted_total)
summary(fitted_total, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary


ds_total <- list(
  prv = prv,
  noti = noti,
  exo = exo
)


save(fitted_total, ds_total, file = "out/Full/KEN/Total.rdata")
# fitted_total, ds_total

#### Marginalised Age ----
prv <- prevalence %>% 
  group_by(Year, Agp) %>%
  summarise(N = round(sum(N)), Asym = round(sum(Asym)), SymSn = round(sum(SymSn)), SymSp = round(sum(SymSp))) %>%
  arrange(Agp)


noti <- notification %>%
  group_by(Year, Agp) %>%
  summarise(n_all = round(sum(n_all)), 
            n_sp = round(sum(n_sp)), 
            n_sn = round(sum(n_sn)), 
            Pop = round(sum(Pop))) %>%
  arrange(Agp)

mor <- mortality %>%
  group_by(Year, Agp) %>%
  summarise(DeaR = weighted.mean(DeaR, Noti)) %>%
  arrange(Year, Agp)


exo <- get_exo_anp(mor, untr_a = "no_untr", untr_s = "full", bg_death = T, pr_fp = 0)


yrs <- sort(unique(noti$Year))
n_t <- length(yrs)

dat <- list(
  YearSurveyed = prv$Year[1],
  N = round(prv$N),
  Asym = round(prv$Asym),
  Sn = round(prv$SymSn),
  Sp = round(prv$SymSp),
  Years = yrs,
  Pop = t(matrix(round(noti$Pop), n_t)),
  NotiSn = t(matrix(round(noti$n_sn), n_t)),
  NotiSp = t(matrix(round(noti$n_sp), n_t)),
  n_t = n_t,
  n_gp = nrow(prv)
)

model <- readRDS(file = "stan/m1_duration_free.rds")
fitted_age <- sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
check_divergences(fitted_age)
summary(fitted_age, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary


model <- readRDS(file = "stan/m1_reg.rds")
fitted_age_cov <- sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
check_divergences(fitted_age_cov)
summary(fitted_age_cov, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr", "lrr_cs_sn", "lrr_cs_sp", "lrr_sym"))$summary


ds_age <- list(
  prv = prv,
  noti = noti,
  exo = exo
)

save(fitted_age, fitted_age_cov, ds_age, file = "out/Full/KEN/Age.rdata")


#### Marginalised Sex ----
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

mor <- mortality %>%
  group_by(Year, Sex) %>%
  summarise(DeaR = weighted.mean(DeaR, Noti)) %>%
  arrange(Year, Sex)

exo <- get_exo_anp(mor, untr_a = "no_untr", untr_s = "full", bg_death = T, pr_fp = 0)


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


model <- readRDS(file = "stan/m1_duration_free.rds")
fitted_sex <- sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
check_divergences(fitted_sex)
summary(fitted_sex, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary


model <- readRDS(file = "stan/m1_reg.rds")
fitted_sex_cov <- sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
check_divergences(fitted_sex_cov)
summary(fitted_sex_cov, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr", "lrr_cs_sn", "lrr_cs_sp", "lrr_sym"))$summary


ds_sex <- list(
  prv = prv,
  noti = noti,
  exo = exo
)

save(fitted_sex, fitted_sex_cov, ds_sex, file = "out/Full/KEN/Sex.rdata")



#### Marginalised HIV ----
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

mor <- mortality %>% 
  left_join(tibble(Year = 2016, HIV = c("HIV", "NonHIV"))) %>%
  mutate(DeaR = ifelse(HIV == "HIV", DeaR + 21/1300, DeaR)) %>%
  group_by(Year, HIV) %>%
  summarise(DeaR = weighted.mean(DeaR, Noti)) %>%
  arrange(Year, HIV)

exo <- get_exo_anp_hiv(mor, untr_a = "no_untr", untr_s = "full", untr_hiv = "as_nonhiv", bg_death = T, pr_fp = 0)


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


model <- readRDS(file = "stan/m1_duration_free.rds")
fitted_hiv <- sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
check_divergences(fitted_hiv)
summary(fitted_hiv, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

model <- readRDS(file = "stan/m1_reg.rds")
fitted_hiv_cov <- sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
check_divergences(fitted_hiv_cov)
summary(fitted_hiv_cov, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr", "lrr_cs_sn", "lrr_cs_sp", "lrr_sym"))$summary


ds_hiv <- list(
  prv = prv,
  noti = noti,
  exo = exo
)

save(fitted_hiv, fitted_hiv_cov, ds_hiv, file = "out/Full/KEN/HIV.rdata")



#### Multivariate ----

prv <- prevalence %>% 
  mutate(N = round(N), Asym = round(Asym), SymSn = round(SymSn), SymSp = round(SymSp)) %>%
  arrange(Sex, Agp, HIV)


noti <- notification %>%
  arrange(Year, Sex, Agp, HIV)

mor <- mortality %>% 
  left_join(tibble(Year = 2016, HIV = c("HIV", "NonHIV"))) %>%
  mutate(DeaR = ifelse(HIV == "HIV", DeaR + 21/1300, DeaR)) %>%
  arrange(Year, Sex, Agp, HIV)

#exo <- get_exo_anp(mor, untr_a = "no_untr", untr_s = "full", bg_death = T, pr_fp = 0)
exo <- get_exo_anp_hiv(mor, untr_a = "no_untr", untr_s = "full", untr_hiv = "as_nonhiv", bg_death = T, pr_fp = 0)


yrs <- sort(unique(noti$Year))
n_t <- length(yrs)

dat <- list(
  YearSurveyed = prv$Year[1],
  N = round(prv$N),
  Asym = round(prv$Asym),
  Sn = round(prv$SymSn),
  Sp = round(prv$SymSp),
  Years = yrs,
  Pop = round(matrix(noti$Pop, nrow(prv), n_t)),
  NotiSn = round(matrix(noti$n_sn, nrow(prv), n_t)),
  NotiSp = round(matrix(noti$n_sp, nrow(prv), n_t)),
  Cov = cbind(prv$Sex == "Male", 
              prv$Agp == "[25,35)", 
              prv$Agp == "[35,45)", prv$Agp == "[45,55)", 
              prv$Agp == "[55,65)", prv$Agp == "[65,Inf)", 
              prv$HIV == "HIV") + 0,
  n_t = n_t,
  n_gp = nrow(prv),
  n_cov = 7
)


model <- readRDS(file = "stan/m2_cov.rds")
fitted_full <- sampling(model, data = c(dat, exo), iter = 3000, warmup = 2000, chain = 3)

summary(fitted_full, pars = c("r_sym", "r_det_sn", "r_det_sp"))$summary
summary(fitted_full, pars = c("lrr_cs_sn", "lrr_cs_sp", "lrr_sym"))$summary

ds_full <- list(
  prv = prv,
  noti = noti,
  exo = exo
)

save(fitted_full, ds_full, file = "out/Full/KEN/Cov.rdata")


save(fitted_age, fitted_age_cov, ds_age,
     fitted_sex, fitted_sex_cov, ds_sex,
     fitted_hiv, fitted_hiv_cov, ds_hiv,
     fitted_full, ds_full,
     fitted_total, 
     file = "out/Full/KEN/All.rdata")

