rm(list = ls())
library(tidyverse)

source("data/country_list.R")

frm_mci <- function(ds, digits = 1) {
  sprintf("%s (%s-%s)", 
          format(mean(ds, na.rm = T), nsmall = digits, digits = digits), 
          format(quantile(ds, 0.025, na.rm = T), nsmall = digits, digits = digits), 
          format(quantile(ds, 0.975, na.rm = T), nsmall = digits, digits = digits))
}

frm_pci <- function(ds, digits = NULL) {
  sprintf("%s (%s-%s)", 
          scales::percent(mean(ds, na.rm = T), accuracy = digits), 
          scales::percent(quantile(ds, 0.025, na.rm = T), accuracy = digits), 
          scales::percent(quantile(ds, 0.975, na.rm = T), accuracy = digits))
}


load("out/summary_data.rdata")
load("out/Durations_all.rdata")


Tab_Dur <- 
  Durations$Durations_All %>%
  group_by(Country, ISO) %>%
  summarise(PN_All = frm_mci(PN_All, 1),
            Delay_All = frm_mci(DelayA * 12, 1),
            CDR_All = frm_pci(CDR_A),
            Inc_All = frm_mci(IncR_All * 1E5, 0),
            
            PN_Sym = frm_mci(PN_S, 1),
            Delay_Sym = frm_mci(DelayS * 12, 1),
            CDR_Sym = frm_pci(CDR_S),
            Inc_Sym = frm_mci(IncR_S * 1E5, 0),
            ) %>%
  left_join(Tab_dat %>% mutate(
    Pop = round(Pop * 1E-6, 1),
    Inc_WHO = sprintf("%s (%s-%s)", round(IncR_M * 1E5), round(IncR_L * 1E5), round(IncR_U * 1E5)), 
    Noti_WHO = round(CNR * 1E5),
    CDR_WHO = sprintf("%s (%s-%s)", 
                      scales::percent(CDR_M / 100, accuracy = 1),
                      scales::percent(CDR_L / 100, accuracy = 1),
                      scales::percent(CDR_U / 100, accuracy = 1))
    ) %>%
    select(Country, Pop, Noti_WHO, CDR_WHO, Inc_WHO)) %>%
  mutate(Country = countries_lab[Country])
  
  
t(Tab_Dur)

write.csv(t(Tab_Dur), file = "docs/tabs/Burden.csv")


