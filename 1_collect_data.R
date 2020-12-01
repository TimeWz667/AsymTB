library(tidyverse)


source("data/country_list.R")



Tab_dat <- bind_rows(lapply(names(countries), function(iso) {
  country <- countries[iso]
  
  iso <- glue::as_glue(iso)
  load("data/Input_" + iso + ".rdata")
  
  
  prevalence %>%
    mutate(Country = country) %>%
    group_by(Country, Year) %>%
    summarise(N = sum(N),
              Asym = sum(Asym),
              Sym = sum(Sym),
              Sn = sum(Sn), 
              Sp = sum(Sp)) %>%
    mutate(Prv = Sn + Sp, PrvPrSp = Sn / Prv, PrvPrAsym = Asym / Prv,
           PrvP_M = Prv / N, PrvP_L = qbinom(0.025, N, Prv / N) / N, PrvP_U = qbinom(0.975, N, Prv / N) / N) %>%
    select(Country, YearSurvey = Year, PrvPrSp, PrvPrAsym, starts_with("PrvP")) %>%
    left_join(notification %>% 
                mutate(Country = country) %>%
                group_by(Country, Year) %>%
                summarise(NotiN = sum(n_all), n_sn = sum(n_sn), n_sp = sum(n_sp), Pop = sum(Pop)) %>%
                filter(Year == 2019) %>%
                mutate(NotiPrSp = n_sp / NotiN)) %>%
    left_join(incidence %>% filter(Sex == "Total") %>%
                rename(Inc_M = m, Inc_L = l, Inc_U = u) %>%
                select(Country, Year, Inc_M, Inc_L, Inc_U)) %>%
    left_join(cdr %>% filter(Year == 2018) %>%
              select(Country, CDR_M = m, CDR_L = l, CDR_U = u)) %>%
    mutate(CNR_M = NotiN / Pop, CNR_L = qbinom(0.025, Pop, NotiN / Pop) / Pop, 
           CNR_U = qbinom(0.975, Pop, NotiN / Pop) / Pop,
           IncR_M = Inc_M / Pop, IncR_L = Inc_L / Pop, IncR_U = Inc_U / Pop,
           CDR_L, CDR_M, CDR_U)
}))


t(Tab_dat)

save(Tab_dat, file = "out/summary_data.rdata")

