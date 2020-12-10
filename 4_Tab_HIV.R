library(tidyverse)
library(rstan)


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


pars <- c("dur_a", "p_sp", "Delay_Sn", "Delay_Sp", "IncN_A", "IncN_S", 
          "CDR_A", "CDR_S", "PrvN_A", "PrvN_Sn", "PrvN_Sp")


load("out/Sens/Post_HIV.rdata")

#### Replace by the main fitting (all data and assumptions are identical)
load("out/Full/KEN/HIV.rdata")
fitted_ken$H2 <- fitted_hiv

load("out/Full/BLT/HIV.rdata")
fitted_blt$H2 <- fitted_hiv


load("data/Input_BLT.rdata")
noti_blt <- notification %>% 
  group_by(Year, HIV) %>% 
  summarise(Pop = sum(Pop)) %>% 
  filter(Year == max(notification$Year))

load("data/Input_KEN_Full.rdata")
noti_ken <- notification %>% 
  group_by(Year, HIV) %>% 
  summarise(Pop = sum(Pop)) %>% 
  filter(Year == max(notification$Year))



hiv_blt <- lapply(fitted_blt, function(fit) {
  with(extract(fit, pars = pars), {
    tibble(
      HIV = rep(c("PLHIV", "Non-HIV"), each = nrow(dur_a)),
      Pop = rep(noti_blt$Pop, each = nrow(dur_a)),
      Dur_A = c(dur_a),
      p_sp = rep(p_sp, 2),
      DelaySn = c(Delay_Sn),
      DelaySp = c(Delay_Sp),
      IncN_All = c(IncN_A),
      IncN_S = c(IncN_S),
      CDR_All = c(CDR_A),
      CDR_S = c(CDR_S),
      PrvN_All = c(PrvN_A + PrvN_Sn + PrvN_Sp),
      PrvN_S = c(PrvN_Sn + PrvN_Sp)
    )
  })
})


hiv_ken <- lapply(fitted_ken, function(fit) {
  with(extract(fit, pars = pars), {
    tibble(
      HIV = rep(c("PLHIV", "Non-HIV"), each = nrow(dur_a)),
      Pop = rep(noti_ken$Pop, each = nrow(dur_a)),
      Dur_A = c(dur_a),
      p_sp = rep(p_sp, 2),
      DelaySn = c(Delay_Sn),
      DelaySp = c(Delay_Sp),
      IncN_All = c(IncN_A),
      IncN_S = c(IncN_S),
      CDR_All = c(CDR_A),
      CDR_S = c(CDR_S),
      PrvN_All = c(PrvN_A + PrvN_Sn + PrvN_Sp),
      PrvN_S = c(PrvN_Sn + PrvN_Sp)
    )
  })
})

for (key in names(hiv_blt)) {
  hiv_blt[[key]] <- hiv_blt[[key]] %>% mutate(Case = key, Location = "Blantyre")
  hiv_ken[[key]] <- hiv_ken[[key]] %>% mutate(Case = key, Location = "Kenya")
}

Tab_HIV <- bind_rows(c(
  hiv_blt,
  hiv_ken
)) %>% 
  mutate(DelaySym = (1 - p_sp) * DelaySn + p_sp * DelaySp,
         DelayAll = Dur_A + DelaySym) %>%
  group_by(Location, Case, HIV) %>%
  summarise(
    DurA = frm_mci(Dur_A * 12),
    DelaySn = frm_mci(DelaySn * 12),
    DelaySp = frm_mci(DelaySp * 12),
    DelaySym = frm_mci(DelaySym * 12),
    DelayAll = frm_mci(DelayAll * 12),
    
    PrvAll = frm_mci(PrvN_All / Pop * 1E5, digits = 0),
    CDR_All = frm_pci(CDR_All, 1),
    IncAll = frm_mci(IncN_All / Pop * 1E5, digits = 0),
    
    PrvSym = frm_mci(PrvN_S / Pop * 1E5, digits = 0),
    CDR_Sym = frm_pci(CDR_S, 1),
    IncSym = frm_mci(IncN_S / Pop * 1E5, digits = 0)
  )



write.csv(t(Tab_HIV %>% arrange(Location, HIV, Case) %>% filter(HIV == "PLHIV" | Case == "H2")), file = "docs/tabs/HIV.csv")





