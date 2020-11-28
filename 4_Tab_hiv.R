rm(list = ls())
library(tidyverse)
library(rstan)

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


pars <- c("r_sym", "r_det_sp", "r_det_sn", 
          "Gap_A", "Gap_S", "Gap_Sp", "Delay_Sn", "Delay_Sp", "CDR_A", "CDR_S", "CDR_Sp",
          "IncN_A", "IncN_S", "NotiN_Sn", "NotiN_Sp", "dur_a", "dur_sn", "dur_sp")
pars <- apply(expand.grid(pars, c("[1]", "[2]")), 1, function(x) paste0(x, collapse = ""))
pars <- c(pars, "p_sp")




Tab_HIV <- bind_rows(list(
  local({
    load("out/Full/BLT/HIV.rdata")
    
    n_iter <- length(extract(fitted_hiv, pars = "r_tr")$r_tr)
    
    pop <- ds_hiv$noti %>% filter(Year == 2018.875) %>% ungroup() %>% 
      select(Pop, HIV) %>% 
      mutate(Location = "Blantyre") 
    
    
    as_tibble(extract(fitted_hiv, pars = pars)) %>%
      mutate(Key = 1:n_iter, "p_sp[1]" = p_sp, "p_sp[2]" = p_sp) %>%
      select(-p_sp) %>%
      pivot_longer(-Key) %>%
      tidyr::extract(name, c("Index", "HIV"), "(\\w+)\\[(\\d)\\]") %>%
      mutate(HIV = ifelse(HIV == 1, "HIV", "NonHIV"), Location = "Blantyre") %>%
      pivot_wider(names_from = Index, values_from = value) %>%
      group_by(Location, HIV) %>%
      left_join(pop) 
    
  }),
  local({
    load("out/Full/KEN/HIV.rdata")
    
    n_iter <- length(extract(fitted_hiv, pars = "r_tr")$r_tr)
    
    pop <- ds_hiv$noti %>% filter(Year == 2019) %>% ungroup() %>% 
      select(Pop, HIV) %>% 
      mutate(Location = "Kenya") 
    
    
    as_tibble(extract(fitted_hiv, pars = pars)) %>%
      mutate(Key = 1:n_iter, "p_sp[1]" = p_sp, "p_sp[2]" = p_sp) %>%
      select(-p_sp) %>%
      pivot_longer(-Key) %>%
      tidyr::extract(name, c("Index", "HIV"), "(\\w+)\\[(\\d)\\]") %>%
      mutate(HIV = ifelse(HIV == 1, "HIV", "NonHIV"), Location = "Kenya") %>%
      pivot_wider(names_from = Index, values_from = value) %>%
      group_by(Location, HIV) %>%
      left_join(pop)
    
  })
)) %>%
  summarise(
    Inc_A = frm_mci(IncN_A / Pop * 1E5, 0),
    #R_Sym = frm_mci(r_sym), 
    Duration_A = frm_mci(dur_a * 12),
    Delay_A = frm_mci((dur_a + p_sp * Delay_Sp + (1 - p_sp) * Delay_Sn) * 12),
    CDR_A = frm_pci(CDR_A), 
    
    Inc_S = frm_mci(IncN_S / Pop * 1E5, 0),
    #R_CS_Sn = frm_mci(r_det_sn), 
    #R_CS_Sp = frm_mci(r_det_sp), 
    Delay_Sn = frm_mci(Delay_Sn * 12), 
    Delay_Sp = frm_mci(Delay_Sp * 12),
    CDR_S = frm_pci(CDR_S), 
  )


t(Tab_HIV)[, 4:1]
write.csv(t(Tab_HIV)[, 4:1], file = "docs/tabs/HIV.csv")

