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


pars <- c("dur_a", "dur_sp", "dur_sn", "Delay_Sn", "Delay_Sp", "IncN_A", "IncN_S", 
          "CDR_A", "CDR_S", "PrvN_A", "PrvN_Sn", "PrvN_Sp")
pars <- apply(expand.grid(pars, c("[1]", "[2]")), 1, function(x) paste0(x, collapse = ""))


order_out <- c("Dur_A", "Dur_All", "CDR_All", "Inc_All", "Dur_Sym", "CDR_Sym", "Inc_Sym")


##################################################################
##### Sensitivity analysis: different definition of symptoms #####
##################################################################

#### Load data
load("out/Sens/Post_Sym.rdata")


Tab_Sym <- bind_rows(lapply(names(fitted_sym), function(key) {
  fi <- fitted_sym[[key]]

  if (startsWith(key, "B")) {
    pop <- population$BLT
  } else {
    pop <- population$KEN
  }
  
  p_sp <- extract(fi, "p_sp")$p_sp
  n_iter <- length(p_sp)


  as_tibble(extract(fi, pars = pars)) %>%
    mutate(Key = 1:n_iter, `p_sp[1]` = p_sp, `p_sp[2]` = p_sp) %>%
    pivot_longer(-Key) %>%
    tidyr::extract(name, c("Index", "Sex"), "(\\w+)\\[(\\d)\\]") %>%
    mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>%
    pivot_wider(names_from = Index, values_from = value) %>%
    left_join(pop) %>%
    rename(DurA = dur_a, DurSymSn = dur_sn) %>%
    mutate(DurSymSp = dur_sp, DurSym = (1 - p_sp) * Delay_Sn + p_sp * Delay_Sp,
           DurAll = DurA + DurSym) %>%
    group_by(Key) %>%
    summarise(Dur_A = weighted.mean(DurA, IncN_A),
              Dur_All = weighted.mean(DurAll, IncN_A),
              CDR_All = weighted.mean(CDR_A, IncN_A),
              Inc_All = sum(IncN_A) / sum(Pop),
              Dur_Sym = weighted.mean(DurSym, IncN_S),
              CDR_Sym = weighted.mean(CDR_S, IncN_S),
              Inc_Sym = sum(IncN_S) / sum(Pop)) %>%
    mutate(Case = key, Location = ifelse(startsWith(key, "B"), "Blantyre", "Kenya"))
})) %>%
  group_by(Location, Case) %>%
  summarise(
    Dur_A = frm_mci(Dur_A * 12, 2),
    Dur_All = frm_mci(Dur_All * 12, 2),
    CDR_All = frm_pci(CDR_All),
    Inc_All = frm_mci(Inc_All * 1E5, 0),
    
    Dur_Sym = frm_mci(Dur_Sym * 12, 2),
    CDR_Sym = frm_pci(CDR_Sym),
    Inc_Sym = frm_mci(Inc_Sym * 1E5, 0)
  ) %>%
  pivot_longer(c(-Location, -Case), names_to = "Variable") %>%
  pivot_wider(Variable, names_from = Case, values_from = value) %>%
  mutate(Variable = factor(Variable, levels = order_out)) %>%
  arrange(Variable)


write.csv(Tab_Sym, file = "docs/tabs/Sym.csv")



#######################################################################
##### Sensitivity analysis: different HIV specific TB death rates #####
#######################################################################

#### Load data
load("out/Sens/Post_HIV.rdata")


population <- list(
  BLT = local({
    load("data/Input_BLT.rdata")
    
    notification %>%
      group_by(Year, HIV) %>%
      summarise(Pop = round(sum(Pop))) %>%
      mutate(HIV = ifelse(HIV == "HIV", "PLHIV", "Non-HIV")) %>%
      arrange(HIV) %>% filter(Year == max(notification$Year))
  }),
  KEN = local({
    load("data/Input_KEN_Full.rdata")
    
    notification %>%
      group_by(Year, HIV) %>%
      summarise(Pop = round(sum(Pop))) %>%
      mutate(HIV = ifelse(HIV == "HIV", "PLHIV", "Non-HIV")) %>%
      arrange(HIV) %>% filter(Year == max(notification$Year))
  })
)


Tab_HIV_BLT <- bind_rows(lapply(names(fitted_blt), function(key) {
  fi <- fitted_blt[[key]]
  
  pop <- population$BLT
  
  p_sp <- extract(fi, "p_sp")$p_sp
  n_iter <- length(p_sp)
  
  
  as_tibble(extract(fi, pars = pars)) %>%
    mutate(Key = 1:n_iter, `p_sp[1]` = p_sp, `p_sp[2]` = p_sp) %>%
    pivot_longer(-Key) %>%
    tidyr::extract(name, c("Index", "HIV"), "(\\w+)\\[(\\d)\\]") %>%
    mutate(HIV = ifelse(HIV == 1, "PLHIV", "Non-HIV")) %>%
    pivot_wider(names_from = Index, values_from = value) %>%
    left_join(pop) %>%
    rename(Dur_A = dur_a, DurSymSn = dur_sn, CDR_Sym = CDR_S, CDR_All = CDR_A) %>%
    mutate(DurSymSp = dur_sp, Dur_Sym = (1 - p_sp) * Delay_Sn + p_sp * Delay_Sp,
           Dur_All = Dur_A + Dur_Sym,
           Inc_Sym = IncN_S / Pop, Inc_All = IncN_A / Pop,
           Case = key, Location = "Blantyre")
}))


Tab_HIV_KEN <- bind_rows(lapply(names(fitted_ken), function(key) {
  fi <- fitted_ken[[key]]
  
  pop <- population$KEN
  
  p_sp <- extract(fi, "p_sp")$p_sp
  n_iter <- length(p_sp)
  
  
  as_tibble(extract(fi, pars = pars)) %>%
    mutate(Key = 1:n_iter, `p_sp[1]` = p_sp, `p_sp[2]` = p_sp) %>%
    pivot_longer(-Key) %>%
    tidyr::extract(name, c("Index", "HIV"), "(\\w+)\\[(\\d)\\]") %>%
    mutate(HIV = ifelse(HIV == 1, "PLHIV", "Non-HIV")) %>%
    pivot_wider(names_from = Index, values_from = value) %>%
    left_join(pop) %>%
    rename(Dur_A = dur_a, DurSymSn = dur_sn, CDR_Sym = CDR_S, CDR_All = CDR_A) %>%
    mutate(DurSymSp = dur_sp, Dur_Sym = (1 - p_sp) * Delay_Sn + p_sp * Delay_Sp,
           Dur_All = Dur_A + Dur_Sym,
           Inc_Sym = IncN_S / Pop, Inc_All = IncN_A / Pop,
           Case = key, Location = "Kenya")
}))


Tab_HIV <- bind_rows(Tab_HIV_BLT, Tab_HIV_KEN) %>%
  group_by(Location, Case, HIV) %>%
  summarise(
    Dur_A = frm_mci(Dur_A * 12, 2),
    Dur_All = frm_mci(Dur_All * 12, 2),
    CDR_All = frm_pci(CDR_All),
    Inc_All = frm_mci(Inc_All * 1E5, 0),
    
    Dur_Sym = frm_mci(Dur_Sym * 12, 2),
    CDR_Sym = frm_pci(CDR_Sym),
    Inc_Sym = frm_mci(Inc_Sym * 1E5, 0)
  ) %>%
  pivot_longer(c(-Location, -Case, -HIV), names_to = "Variable") %>%
  filter(Case == "H1" | HIV == "PLHIV") %>%
  mutate(name = paste(Location, Case, HIV, sep=":")) %>%
  pivot_wider(Variable, names_from = name, values_from = value) %>%
  mutate(Variable = factor(Variable, levels = order_out)) %>%
  arrange(Variable)

write.csv(Tab_HIV, file = "docs/tabs/HIV.csv")
