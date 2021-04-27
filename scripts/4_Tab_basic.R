library(tidyverse)

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


load("out/summary_data.rdata")
load("out/Durations_All.rdata")


Tab_Dur <- Durations$Durations_All %>%
  group_by(Country, ISO) %>%
  select(-CNR) %>%
  left_join(Tab_dat %>% select(Country, CNR = CNR_M)) %>%
  summarise(#Prv_All = frm_mci(PrvN_All / Pop * 1E5, 0),
            PN_All = frm_mci(PrvN_All / Pop / CNR * 12, 1),
            Delay_All = frm_mci(Delay_All * 12, 1),
            CDR_All = frm_pci(CDR_A),
            Inc_All = frm_mci(IncR_All * 1E5, 0),
            #Prv_Sym = frm_mci((PrvN_All - PrvN_A) / Pop * 1E5, 0),
            PN_Sym = frm_mci((PrvN_All - PrvN_A) / Pop / CNR * 12, 1),
            Delay_SymCs = frm_mci(DurS * 12, 1),
            Delay_Sym = frm_mci(Delay_S * 12, 1),
            CDR_Sym = frm_pci(CDR_S),
            Inc_Sym = frm_mci(IncR_S * 1E5, 0),
            ADR = frm_pci(ADR, 0.1),
            PropDurA = frm_pci(DurA / TTN),
            PropDurSymS = frm_pci(DurS / (DurS + DurC)),
            PrSC = frm_pci(PrSC),
            PrMor = frm_pci(PrMor)
  ) %>%
  left_join(Tab_dat %>% mutate(
    PopSize = sprintf("%.1f", Pop * 1E-6),
    YearSurvey = YearSurvey,
    PrvP_Survey = sprintf("%.0f (%.0f \u2013 %.0f)", round(PrvP_M * 1E5), round(PrvP_L * 1E5), round(PrvP_U * 1E5)),
    Noti_WHO = sprintf("%.0f (%.0f \u2013 %.0f)", round(CNR_M * 1E5), round(CNR_L * 1E5), round(CNR_U * 1E5)),
    CDR_WHO = sprintf("%.0f (%.1f \u2013 %.0f)", CDR_M, CDR_L, CDR_U),
    Inc_WHO = sprintf("%.0f (%.0f \u2013 %.0f)", round(IncR_M * 1E5), round(IncR_L * 1E5), round(IncR_U * 1E5))
  ) %>%
    select(Country, PopSize, YearSurvey, PrvP_Survey, Noti_WHO, CDR_WHO, Inc_WHO)) %>%
  mutate(Country = countries_lab[Country],
         Delay_SymCs = ifelse(Delay_SymCs == Delay_Sym, "N.A.", Delay_SymCs)) %>%
  relocate(PopSize, YearSurvey, PrvP_Survey, Noti_WHO, CDR_WHO, Inc_WHO ,.after = ISO)
  
  
t(Tab_Dur)

write.csv(t(Tab_Dur), file = "docs/tabs/Burden.csv")


