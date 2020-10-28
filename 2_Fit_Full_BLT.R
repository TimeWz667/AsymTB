library(tidyverse)
library(rstan)

options(mc.cores = min(5, parallel::detectCores()))



source("R/get_exo.R")


load("data/Input_Blantyre.rdata")


#### Aggregated ----
prv <- prevalence
noti <- notification %>%
  group_by(Year) %>%
  summarise(n_all = sum(n_all), n_sp = sum(n_sp), n_sn = sum(n_sn), Pop = sum(Pop))

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

dr <- mortality %>%
  group_by(Year) %>%
  summarise(DeaR = sum(DeaR * Pop) / sum(Pop))
dr <- dr$DeaR


model <- readRDS(file = "stan/m0.rds")

exo <- get_exo(F, dr, untr = "no_untr", bg_death = T)
fitted1 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted1, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

exo <- get_exo(F, dr, untr = "as_sn", bg_death = T)
fitted2 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted2, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

exo <- get_exo(F, dr, untr = "full", bg_death = T)
fitted3 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted3, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

dataset <- list(
  prv = prv,
  noti = noti,
  exo = exo
)


save(fitted1, fitted2, fitted3, dataset, file = "out/Full/BLT/Total.rdata")


#### Marginalised Age ----
prv <- prevalence %>% 
  group_by(Year, Agp) %>%
  summarise(N = round(sum(N)), Asym = round(sum(Asym)), SymSn = round(sum(SymSn)), SymSp = round(sum(SymSp))) %>%
  arrange(Agp)
  
  
noti <- notification %>%
  group_by(Year, Agp) %>%
  summarise(n_all = sum(n_all), n_sp = sum(n_sp), n_sn = sum(n_sn), Pop = sum(Pop)) %>%
  arrange(Agp)

dr <- mortality %>%
  group_by(Year, Agp) %>%
  summarise(DeaR = sum(DeaR * Pop) / sum(Pop)) %>%
  arrange(Agp)
dr <- dr$DeaR


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


exo <- get_exo(rep(F, nrow(prv)), dr, untr = "full", bg_death = T)


model <- readRDS(file = "stan/m1_uni.rds")
fitted0 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted0, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

model <- readRDS(file = "stan/m1_duration_free.rds")
fitted1 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted1, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

model <- readRDS(file = "stan/m1_fixed.rds")
fitted2 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted2, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

model <- readRDS(file = "stan/m1_reg.rds")
fitted3 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted3, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr", "lrr_cs_sn", "lrr_cs_sp", "lrr_sym"))$summary

wts <- loo::loo_model_weights(list(loo(fitted0, "log_lik_pr"), 
                                   loo(fitted1, "log_lik_pr"), 
                                   loo(fitted2, "log_lik_pr"),
                                   loo(fitted3, "log_lik_pr")))

dataset <- list(
  prv = prv,
  noti = noti,
  exo = exo
)

save(fitted0, fitted1, fitted2, fitted3, dataset, file = "out/Full/BLT/Age.rdata")


#### Marginalised Sex ----
prv <- prevalence %>% 
  group_by(Year, Sex) %>%
  summarise(N = round(sum(N)), Asym = round(sum(Asym)), SymSn = round(sum(SymSn)), SymSp = round(sum(SymSp))) %>%
  arrange(Sex)


noti <- notification %>%
  group_by(Year, Sex) %>%
  summarise(n_all = sum(n_all), n_sp = sum(n_sp), n_sn = sum(n_sn), Pop = sum(Pop)) %>%
  arrange(Sex)

dr <- mortality %>%
  group_by(Year, Sex) %>%
  summarise(DeaR = sum(DeaR * Pop) / sum(Pop)) %>%
  arrange(Sex)
dr <- dr$DeaR


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


exo <- get_exo(rep(F, nrow(prv)), dr, untr = "full", bg_death = T)


model <- readRDS(file = "stan/m1_uni.rds")
fitted0 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted0, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

model <- readRDS(file = "stan/m1_duration_free.rds")
fitted1 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted1, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

model <- readRDS(file = "stan/m1_fixed.rds")
fitted2 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted2, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

model <- readRDS(file = "stan/m1_reg.rds")
fitted3 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted3, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr", "lrr_cs_sn", "lrr_cs_sp", "lrr_sym"))$summary


wts <- loo::loo_model_weights(list(loo(fitted0, "log_lik_pr"), 
                                   loo(fitted1, "log_lik_pr"), 
                                   loo(fitted2, "log_lik_pr"),
                                   loo(fitted3, "log_lik_pr")))

dataset <- list(
  prv = prv,
  noti = noti,
  exo = exo
)

save(fitted0, fitted1, fitted2, fitted3, dataset, file = "out/Full/BLT/Sex.rdata")



#### Marginalised HIV ----
prv <- prevalence %>% 
  group_by(Year, HIV) %>%
  summarise(N = round(sum(N)), Asym = round(sum(Asym)), SymSn = round(sum(SymSn)), SymSp = round(sum(SymSp))) %>%
  arrange(HIV)


noti <- notification %>%
  group_by(Year, HIV) %>%
  summarise(n_all = sum(n_all), n_sp = sum(n_sp), n_sn = sum(n_sn), Pop = sum(Pop)) %>%
  arrange(HIV)

dr <- mortality %>%
  group_by(Year) %>%
  summarise(DeaR = sum(DeaR * Pop) / sum(Pop))
dr <- dr$DeaR


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


exo <- get_exo(prv$HIV, rep(dr, nrow(prv)), untr = "full", bg_death = T)


model <- readRDS(file = "stan/m1_uni.rds")
fitted0 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted0, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

model <- readRDS(file = "stan/m1_duration_free.rds")
fitted1 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted1, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

model <- readRDS(file = "stan/m1_fixed.rds")
fitted2 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted2, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr"))$summary

model <- readRDS(file = "stan/m1_reg.rds")
fitted3 <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)
summary(fitted3, pars = c("r_sym", "r_det_sn", "r_det_sp", "p_sp", "r_tr", "lrr_cs_sn", "lrr_cs_sp", "lrr_sym"))$summary


wts <- loo::loo_model_weights(list(loo(fitted0, "log_lik_noti"), 
                                   loo(fitted1, "log_lik_noti"), 
                                   loo(fitted2, "log_lik_noti"),
                                   loo(fitted3, "log_lik_noti")))

dataset <- list(
  prv = prv,
  noti = noti,
  exo = exo
)

save(fitted0, fitted1, fitted2, fitted3, dataset, file = "out/Full/BLT/HIV.rdata")



#### Multivariate ----

prv <- prevalence %>% 
  mutate(N = round(N), Asym = round(Asym), SymSn = round(SymSn), SymSp = round(SymSp)) %>%
  arrange(Sex, Agp, HIV)


noti <- notification %>%
  arrange(Year, Sex, Agp, HIV)

dr <- mortality %>% arrange(Sex, Agp)
dr <-rep(dr$DeaR, each = 2)


yrs <- sort(unique(noti$Year))
n_t <- length(yrs)

dat <- list(
  YearSurveyed = prv$Year[1],
  N = round(prv$N),
  Asym = round(prv$Asym),
  Sn = round(prv$SymSn),
  Sp = round(prv$SymSp),
  Years = yrs,
  Pop = matrix(noti$Pop, nrow(prv), n_t),
  NotiSn = matrix(noti$n_sn, nrow(prv), n_t),
  NotiSp = matrix(noti$n_sp, nrow(prv), n_t),
  Cov = cbind(prv$Sex == "Male", prv$Agp == "[50,Inf)", prv$HIV == "HIV") + 0,
  n_t = n_t,
  n_gp = nrow(prv),
  n_cov = 3
)


exo <- get_exo(prv$HIV, dr, untr = "full", bg_death = T)

model <- readRDS(file = "stan/m2_cov.rds")
fitted_full <- sampling(model, data = c(dat, exo), iter = 3000, chain = 2)

summary(fitted_full, pars = c("lrr_cs_sn", "lrr_cs_sp", "lrr_sym"))$summary

dataset <- list(
  prv = prv,
  noti = noti,
  exo = exo
)

save(fitted_full, dataset, file = "out/Full/BLT/Cov.rdata")
