library(tidyverse)
library(rstan)

options(mc.cores = min(5, parallel::detectCores()))


source("R/var_exo_anp.R")

###
n_iter = 3000
n_collect = 1000
n_chain = 3


model <- readRDS(file = "stan/m1_duration_free.rds")

### Blantyre ---

fitted_blt <- local({
  load("data/Input_BLT_Sens.rdata")

  noti <- notification %>%
    group_by(Year, Sex) %>%
    summarise(n_all = round(sum(n_all)),
              n_sp = round(sum(n_sp)),
              n_sn = round(sum(n_sn)),
              Pop = round(sum(Pop))) %>%
    arrange(Sex) %>% 
    filter(Year > 2014)
  
  mor <- mortality %>%
    group_by(Year, Sex) %>%
    summarise(DeaR = weighted.mean(DeaR, Noti)) %>%
    arrange(Year, Sex) %>%
    left_join(noti %>%
                group_by(Sex) %>%
                summarise(pr_sp = sum(n_sp) / sum(n_all)))
  
  exo <- get_exo_anp(mor, untr_a = "no_untr", untr_s = "full", bg_death = T, pr_fp = 0)

  yrs <- sort(unique(noti$Year))
  n_t <- length(yrs)
  
  res <- lapply(prevalence_sens, function(prevalence) {
    prv <- prevalence %>% 
      group_by(Year, Sex) %>%
      summarise(N = round(sum(N)), Asym = round(sum(Asym)), SymSn = round(sum(SymSn)), SymSp = round(sum(SymSp))) %>%
      arrange(Sex)
    
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
    
    sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  })

  names(res) <- c("B1", "B2")
  
  res
})


### Kenya ---

fitted_ken <- local({
  load("data/Input_KEN_Sens.rdata")
  
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
    arrange(Year, Sex) %>%
    left_join(noti %>%
                group_by(Sex) %>%
                summarise(pr_sp = sum(n_sp) / sum(n_all)))
  
  exo <- get_exo_anp(mor, untr_a = "no_untr", untr_s = "full", bg_death = T, pr_fp = 0)
  
  yrs <- sort(unique(noti$Year))
  n_t <- length(yrs)
  
  res <- lapply(prevalence_sens[3:5], function(prevalence) {
    prv <- prevalence %>% 
      group_by(Year, Sex) %>%
      summarise(N = round(sum(N)), Asym = round(sum(Asym)), SymSn = round(sum(SymSn)), SymSp = round(sum(SymSp))) %>%
      arrange(Sex)
    
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
    
    sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  })
  
  names(res) <- c("K1", "K3", "K2")
  
  res
})

population <- list(
  BLT = local({
    load("data/Input_BLT_Sens.rdata")
    
    notification %>%
      group_by(Year, Sex) %>%
      summarise(Pop = round(sum(Pop))) %>%
      arrange(Sex) %>% filter(Year == max(notification$Year))
  }),
  KEN = local({
    load("data/Input_KEN_Sens.rdata")
    
    notification %>%
      group_by(Year, Sex) %>%
      summarise(Pop = round(sum(Pop))) %>%
      arrange(Sex) %>% filter(Year == max(notification$Year))
  })
)
### Output

fitted_sym <- c(fitted_blt, fitted_ken)

save(fitted_sym, population, file = "out/Sens/Post_Sym.rdata")

