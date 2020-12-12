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
  
  exos <- list(
    M000 = get_exo_asc(mor, untr_a = "no_untr", untr_s = "no_untr", bg_death = F, pr_fp = 0),
    M100 = get_exo_asc(mor, untr_a = "no_untr", untr_s = "no_untr", bg_death = T, pr_fp = 0),
    M101 = get_exo_asc(mor, untr_a = "no_untr", untr_s = "as_sn", bg_death = T, pr_fp = 0),
    M102 = get_exo_asc(mor, untr_a = "no_untr", untr_s = "avg", bg_death = T, pr_fp = 0),
    M112 = get_exo_asc(mor, untr_a = "as_sn", untr_s = "avg", bg_death = T, pr_fp = 0)
  ) 
  
  model <- readRDS(file = "stan/m3_as_duration_free.rds")

  fitted_mor <- lapply(exos, function(exo) {
    sampling(model, data = c(dat, exo), iter = n_iter, warmup = n_iter - n_collect, chain = n_chain)
  })
  
  save(fitted_mor, dat, file = "out/Sens/Post_Mor_" + iso + ".rdata")
}

