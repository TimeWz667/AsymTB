library(tidyverse)
library(rstan)


source("data/country_list.R")

frm_mci <- function(ds, digits = 1) {
  sprintf("%s (%s \u2013 %s)", 
          format(mean(ds, na.rm = T), nsmall = digits, digits = digits, big.mark = ","), 
          format(quantile(ds, 0.025, na.rm = T), nsmall = digits, digits = digits, big.mark = ","), 
          format(quantile(ds, 0.975, na.rm = T), nsmall = digits, digits = digits, big.mark = ","))
}

frm_pci <- function(ds, digits = NULL) {
  sprintf("%s (%s \u2013 %s)", 
          scales::percent(mean(ds, na.rm = T), accuracy = digits), 
          scales::percent(quantile(ds, 0.025, na.rm = T), accuracy = digits), 
          scales::percent(quantile(ds, 0.975, na.rm = T), accuracy = digits))
}




pars <- c("dur_a", "dur_s", "CDR_A", "IncN_A", "IncN_S")
pars <- apply(expand.grid(pars, c("[1]", "[2]")), 1, function(x) paste0(x, collapse = ""))


Tab.Mor <- bind_rows(lapply(names(countries), function(iso) {
  iso <- glue::as_glue(iso)
  country <- countries[iso]
  
  load("data/Input_" + iso + ".rdata")
  load("out/Sens/Post_Mor_" + iso + ".rdata")
  
  pop <- notification %>%
    filter(Year == 2019) %>%
    group_by(Sex) %>%
    summarise(Pop = sum(Pop)) %>%
    ungroup() %>%
    select(Sex, Pop)
  
  
  tab <- bind_rows(lapply(names(fitted_mor), function(key) {
    fi <- fitted_mor[[key]]
    n_iter <- length(extract(fi, "r_sc")$r_sc)
    
    as_tibble(extract(fi, pars = pars)) %>%
      mutate(Key = 1:n_iter) %>%
      pivot_longer(-Key) %>%
      tidyr::extract(name, c("Index", "Sex"), "(\\w+)\\[(\\d)\\]") %>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>%
      pivot_wider(names_from = Index, values_from = value) %>%
      rename(DurA = dur_a, DurS = dur_s) %>%
      left_join(pop) %>%
      group_by(Key) %>%
      summarise(DurA = weighted.mean(DurA, IncN_A),
                DurS = weighted.mean(DurS, IncN_A),
                CDR = weighted.mean(CDR_A, IncN_A),
                Inc = sum(IncN_A) / sum(Pop)) %>%
      mutate(Case = key, Country = countries_lab[country], ISO = iso)
  })) %>%
    group_by(Case, Country, ISO) %>%
    summarise(
      DurA = frm_mci(DurA * 12, 2),
      DurS = frm_mci(DurS * 12, 2),
      CDR = frm_pci(CDR),
      Inc = frm_mci(Inc * 1E5, 0)
    )
}))




Tab.Mor <- Tab.Mor %>%
  pivot_longer(c(DurA, DurS, CDR, Inc), names_to = "Variable") %>%
  pivot_wider(c(Country, ISO, Variable), names_from = Case, values_from = value) %>%
  mutate(Variable = factor(Variable, levels = c("DurA", "DurS", "CDR", "Inc"))) %>%
  select(-ISO) %>%
  arrange(Variable, Country)


write.csv(Tab.Mor, file = "docs/tabs/Mor.csv")






