library(tidyverse)
library(rstan)

options(mc.cores = min(4, parallel::detectCores()))


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
  MWI = "Malawi", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda", #
  ZMB = "Zambia"
)




i = 4

iso <- names(countries)[i]
country <- countries[i]


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


model <- readRDS(file = "stan/m3_as_fixed.rds")


res_00 <- sampling(model, data = c(dat, get_exo_sym(dr, "no_untr", F)), iter = 3000, chain = 2)
res_10 <- sampling(model, data = c(dat, get_exo_sym(dr, "as_sn", F)), iter = 3000, chain = 2)
res_20 <- sampling(model, data = c(dat, get_exo_sym(dr, "avg", F)), iter = 3000, chain = 2)
res_01 <- sampling(model, data = c(dat, get_exo_sym(dr, "no_untr", T)), iter = 3000, chain = 2)
res_11 <- sampling(model, data = c(dat, get_exo_sym(dr, "as_sn", T)), iter = 3000, chain = 2)
res_21 <- sampling(model, data = c(dat, get_exo_sym(dr, "avg", T)), iter = 3000, chain = 2)





summary(res_00, pars = c("r_sym", "r_det", "adr", "dur_a", "dur_s"))$summary


summary(res_00, pars = c("inc_a[1,9]", "inc_a[2,9]"))$summary[, c(1, 4, 8)]
summary(res_10, pars = c("inc_a[1,9]", "inc_a[2,9]"))$summary[, c(1, 4, 8)]
summary(res_20, pars = c("inc_a[1,9]", "inc_a[2,9]"))$summary[, c(1, 4, 8)]
summary(res_01, pars = c("inc_a[1,9]", "inc_a[2,9]"))$summary[, c(1, 4, 8)]
summary(res_11, pars = c("inc_a[1,9]", "inc_a[2,9]"))$summary[, c(1, 4, 8)]
summary(res_21, pars = c("inc_a[1,9]", "inc_a[2,9]"))$summary[, c(1, 4, 8)]


rr = rbind(
  summary(res_00, pars = c("inc_s[1,9]", "inc_s[2,9]"))$summary[, c(1, 4, 8)],
  summary(res_10, pars = c("inc_s[1,9]", "inc_s[2,9]"))$summary[, c(1, 4, 8)],
  summary(res_20, pars = c("inc_s[1,9]", "inc_s[2,9]"))$summary[, c(1, 4, 8)],
  summary(res_01, pars = c("inc_s[1,9]", "inc_s[2,9]"))$summary[, c(1, 4, 8)],
  summary(res_11, pars = c("inc_s[1,9]", "inc_s[2,9]"))$summary[, c(1, 4, 8)],
  summary(res_21, pars = c("inc_s[1,9]", "inc_s[2,9]"))$summary[, c(1, 4, 8)]
)


colnames(rr) <- c("M", "L", "U")

data.frame(
  Untr = rep(c("A", "B", "C"), each = 2),
  Bg = rep(c("No", "Yes"), each = 6),
  rr) %>%
  ggplot() +
  geom_pointrange(aes(x = Untr, y = M, ymin = L, ymax = U, colour = Bg))




