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



fetch_dur_a <- function(sub) {
  sub <- glue::as_glue(sub)
  
  load(file = "out/Full/" + sub +"/All.rdata")
  
  ext <- extract(fitted_sex, pars = c("r_sym", "dur_a"))
  exo <- ds_sex$exo
  

  tab_sex <- data.table::data.table(
    Variable = "Sex",
    Value = ds_sex$prv$Sex, 
    Duration = apply(ext$dur_a * 12, 2, frm_mci),
    Rate = apply(ext$r_sym, 2, frm_mci, digit = 2),
    RR_uni = c("reference", frm_mci(exp(extract(fitted_sex_cov, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 1]), digit = 2))
  )
  tab_sex
  

  ext <- extract(fitted_age, pars = c("r_sym", "dur_a"))
  exo <- ds_age$exo
  
  n_agp <- length(unique(ds_age$prv$Agp))
  if (n_agp > 2) {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = ds_age$prv$Agp, 
      Duration = apply(ext$dur_a * 12, 2, frm_mci),
      Rate = apply(ext$r_sym, 2, frm_mci, digit = 2),
      RR_uni = c("reference", apply(exp(extract(fitted_age_cov, pars = c("lrr_sym"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2)),
      RR_multi = c("reference", apply(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2))
    )
  } else {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = ds_age$prv$Agp, 
      Duration = apply(ext$dur_a * 12, 2, frm_mci),
      Rate = apply(ext$r_sym, 2, frm_mci, digit = 2),
      RR_uni = c("reference", frm_mci(exp(extract(fitted_age_cov, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
      RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 2]), digit = 2))
    )
  }

  ext <- extract(fitted_hiv, pars = c("r_sym", "dur_a"))
  exo <- ds_hiv$exo
  
  tab_hiv <- data.table::data.table(
    Variable = "HIV",
    Value = rev(ds_hiv$prv$HIV), 
    Duration = rev(apply(ext$dur_a * 12, 2, frm_mci)),
    Rate = rev(apply(ext$r_sym, 2, frm_mci, digit = 2)),
    RR_uni = c("reference", frm_mci(1/exp(extract(fitted_hiv_cov, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, n_agp + 1]), digit = 2))
  )
  tab_hiv
  

  ext <- extract(fitted_total, pars = c("r_sym", "dur_a"))
  
  
  tab_overall <- data.table::data.table(
    Variable = "Overall",
    Value = "", 
    Duration = frm_mci(ext$dur_a * 12),
    Rate = frm_mci(ext$r_sym, digit = 2),
    RR_uni = "",
    RR_multi = ""
  )
  
  rbind(tab_overall, tab_sex, tab_age, tab_hiv)
  
}


fetch_dur_sn <- function(sub) {
  sub <- glue::as_glue(sub)
  
  load(file = "out/Full/" + sub +"/All.rdata")
  
  ext <- extract(fitted_sex, pars = c("r_det_sn", "dur_sn"))
  exo <- ds_sex$exo
  
  
  tab_sex <- data.table::data.table(
    Variable = "Sex",
    Value = ds_sex$prv$Sex, 
    Duration = apply(ext$dur_sn * 12, 2, frm_mci),
    Rate = apply(ext$r_det_sn, 2, frm_mci, digit = 2),
    RR_uni = c("reference", frm_mci(exp(extract(fitted_sex_cov, pars = c("lrr_cs_sn"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sn"))[[1]][, 1]), digit = 2))
  )
  tab_sex
  
  
  ext <- extract(fitted_age, pars = c("r_det_sn", "dur_sn"))
  exo <- ds_age$exo
  
  n_agp <- length(unique(ds_age$prv$Agp))
  if (n_agp > 2) {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = ds_age$prv$Agp, 
      Duration = apply(ext$dur_sn * 12, 2, frm_mci),
      Rate = apply(ext$r_det_sn, 2, frm_mci, digit = 2),
      RR_uni = c("reference", apply(exp(extract(fitted_age_cov, pars = c("lrr_cs_sn"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2)),
      RR_multi = c("reference", apply(exp(extract(fitted_full, c("lrr_cs_sn"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2))
    )
  } else {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = ds_age$prv$Agp, 
      Duration = apply(ext$dur_sn * 12, 2, frm_mci),
      Rate = apply(ext$r_det_sn, 2, frm_mci, digit = 2),
      RR_uni = c("reference", frm_mci(exp(extract(fitted_age_cov, pars = c("lrr_cs_sn"))[[1]][, 2]), digit = 2)),
      RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sn"))[[1]][, 2]), digit = 2))
    )
  }
  
  ext <- extract(fitted_hiv, pars = c("r_det_sn", "dur_sn"))
  exo <- ds_hiv$exo
  
  tab_hiv <- data.table::data.table(
    Variable = "HIV",
    Value = rev(ds_hiv$prv$HIV), 
    Duration = rev(apply(ext$dur_sn * 12, 2, frm_mci)),
    Rate = rev(apply(ext$r_det_sn, 2, frm_mci, digit = 2)),
    RR_uni = c("reference", frm_mci(1/exp(extract(fitted_hiv_cov, pars = c("lrr_cs_sn"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sn"))[[1]][, n_agp + 1]), digit = 2))
  )
  tab_hiv
  
  
  ext <- extract(fitted_total, pars = c("r_det_sn", "dur_sn"))
  
  
  tab_overall <- data.table::data.table(
    Variable = "Overall",
    Value = "", 
    Duration = frm_mci(ext$dur_sn * 12),
    Rate = frm_mci(ext$r_det_sn, digit = 2),
    RR_uni = "",
    RR_multi = ""
  )
  
  rbind(tab_overall, tab_sex, tab_age, tab_hiv)
}


fetch_dur_sp <- function(sub) {
  sub <- glue::as_glue(sub)
  
  load(file = "out/Full/" + sub +"/All.rdata")
  
  ext <- extract(fitted_sex, pars = c("r_det_sp", "dur_sp"))
  exo <- ds_sex$exo
  
  
  tab_sex <- data.table::data.table(
    Variable = "Sex",
    Value = ds_sex$prv$Sex, 
    Duration = apply(ext$dur_sp * 12, 2, frm_mci),
    Rate = apply(ext$r_det_sp, 2, frm_mci, digit = 2),
    RR_uni = c("reference", frm_mci(exp(extract(fitted_sex_cov, pars = c("lrr_cs_sp"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sp"))[[1]][, 1]), digit = 2))
  )
  tab_sex
  
  
  ext <- extract(fitted_age, pars = c("r_det_sp", "dur_sp"))
  exo <- ds_age$exo
  
  n_agp <- length(unique(ds_age$prv$Agp))
  if (n_agp > 2) {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = ds_age$prv$Agp, 
      Duration = apply(ext$dur_sp * 12, 2, frm_mci),
      Rate = apply(ext$r_det_sp, 2, frm_mci, digit = 2),
      RR_uni = c("reference", apply(exp(extract(fitted_age_cov, pars = c("lrr_cs_sp"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2)),
      RR_multi = c("reference", apply(exp(extract(fitted_full, c("lrr_cs_sp"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2))
    )
  } else {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = ds_age$prv$Agp, 
      Duration = apply(ext$dur_sp * 12, 2, frm_mci),
      Rate = apply(ext$r_det_sp, 2, frm_mci, digit = 2),
      RR_uni = c("reference", frm_mci(exp(extract(fitted_age_cov, pars = c("lrr_cs_sp"))[[1]][, 2]), digit = 2)),
      RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sp"))[[1]][, 2]), digit = 2))
    )
  }
  
  ext <- extract(fitted_hiv, pars = c("r_det_sp", "dur_sp"))
  exo <- ds_hiv$exo
  
  tab_hiv <- data.table::data.table(
    Variable = "HIV",
    Value = rev(ds_hiv$prv$HIV), 
    Duration = rev(apply(ext$dur_sp * 12, 2, frm_mci)),
    Rate = rev(apply(ext$r_det_sp, 2, frm_mci, digit = 2)),
    RR_uni = c("reference", frm_mci(1/exp(extract(fitted_hiv_cov, pars = c("lrr_cs_sp"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sp"))[[1]][, n_agp + 1]), digit = 2))
  )
  tab_hiv
  
  
  ext <- extract(fitted_total, pars = c("r_det_sp", "dur_sp"))
  
  
  tab_overall <- data.table::data.table(
    Variable = "Overall",
    Value = "", 
    Duration = frm_mci(ext$dur_sp * 12),
    Rate = frm_mci(ext$r_det_sp, digit = 2),
    RR_uni = "",
    RR_multi = ""
  )
  
  rbind(tab_overall, tab_sex, tab_age, tab_hiv)
}


fetch_dur_total <- function(sub) {
  sub <- glue::as_glue(sub)
  
  load(file = "out/Full/" + sub +"/All.rdata")
  load(file = "out/Full/" + sub +"/Total.rdata")

  dur_total <- with(extract(fitted_total, pars = c("dur_a", "p_sp", "Delay_Sn", "Delay_Sp", "IncN_A", "IncN_S", 
                                                   "CDR_A", "CDR_S", "PrvN_A", "PrvN_Sn", "PrvN_Sp")), {
    tibble(
      Group = "Overall_Overall",
      Pop = (ds_total$noti %>% filter(Year == max(ds_total$noti$Year)))$Pop,
      Dur_A = c(dur_a),
      p_sp = c(p_sp),
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
  
  dur_sex <- with(extract(fitted_sex, pars = c("dur_a", "p_sp", "Delay_Sn", "Delay_Sp", "IncN_A", "IncN_S", 
                                               "CDR_A", "CDR_S", "PrvN_A", "PrvN_Sn", "PrvN_Sp")), {
    tibble(
      Group = rep(paste0("Sex_", ds_sex$prv$Sex), each = nrow(dur_a)),
      Pop = rep((ds_sex$noti %>% filter(Year == max(ds_sex$noti$Year)))$Pop, each = nrow(dur_a)),
      p_sp = rep(p_sp, 2),
      Dur_A = c(dur_a),
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
  
  dur_age <- with(extract(fitted_age, pars = c("dur_a", "p_sp", "Delay_Sn", "Delay_Sp", "IncN_A", "IncN_S", 
                                               "CDR_A", "CDR_S", "PrvN_A", "PrvN_Sn", "PrvN_Sp")), {
    tibble(
      Group = rep(paste0("Age_", ds_age$prv$Agp), each = nrow(dur_a)),
      Pop = rep((ds_age$noti %>% filter(Year == max(ds_age$noti$Year)))$Pop, each = nrow(dur_a)),
      p_sp = rep(p_sp, nrow(ds_age$prv)),
      Dur_A = c(dur_a),
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
  
  dur_hiv <- with(extract(fitted_hiv, pars = c("dur_a", "p_sp", "Delay_Sn", "Delay_Sp", "IncN_A", "IncN_S", 
                                               "CDR_A", "CDR_S", "PrvN_A", "PrvN_Sn", "PrvN_Sp")), {
    tibble(
      Group = rep(paste0("HIV_", ds_hiv$prv$HIV), each = nrow(dur_a)),
      Pop = rep((ds_hiv$noti %>% filter(Year == max(ds_hiv$noti$Year)))$Pop, each = nrow(dur_a)),
      p_sp = rep(p_sp, nrow(ds_hiv$prv)),
      Dur_A = c(dur_a),
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
  
  bind_rows(dur_total, dur_sex, dur_age, dur_hiv) %>%
    mutate(DelaySym = (1 - p_sp) * DelaySn + p_sp * DelaySp,
           DelayAll = Dur_A + DelaySym) %>%
    group_by(Group) %>%
    summarise(
      DurA = frm_mci(Dur_A * 12),
      DelaySn = frm_mci(DelaySn * 12),
      DelaySp = frm_mci(DelaySp * 12),
      DelaySym = frm_mci(DelaySym * 12),
      DelayAll = frm_mci(DelayAll * 12),
      PrvSym = frm_mci(PrvN_S / Pop * 1E5, digits = 0),
      PrvAll = frm_mci(PrvN_All / Pop * 1E5, digits = 0),
      IncSym = frm_mci(IncN_S / Pop * 1E5, digits = 0),
      IncAll = frm_mci(IncN_All / Pop * 1E5, digits = 0),
      CDR_Sym = frm_pci(CDR_S, 1),
      CDR_All = frm_pci(CDR_All, 1),
    ) %>%
    separate(Group, c("Variable", "Value"), "_") %>%
    mutate(Variable = factor(Variable, c("Overall", "Sex", "Age", "HIV")),
           Value = ifelse(Value == "Overall", "", Value),
           Value = factor(Value, c("", "Female", "Male", as.character(ds_age$prv$Agp), "NonHIV", "HIV"))) %>%
    arrange(Variable)
  
}


tab_asym_KEN <- fetch_dur_a("KEN")
tab_sn_KEN <- fetch_dur_sn("KEN")
tab_sp_KEN <- fetch_dur_sp("KEN")
tab_dur_KEN <- fetch_dur_total("KEN")

tab_asym_BLT <- fetch_dur_a("BLT")
tab_sn_BLT <- fetch_dur_sn("BLT")
tab_sp_BLT <- fetch_dur_sp("BLT")
tab_dur_BLT <- fetch_dur_total("BLT")



#### Output ----
write.csv(tab_asym_BLT, file = "docs/tabs/DurAsym_BLT.csv")
write.csv(tab_asym_KEN, file = "docs/tabs/DurAsym_KEN.csv")

write.csv(bind_rows(cbind(Location = "Blantyre", tab_asym_BLT %>% 
                            left_join(tab_dur_BLT %>% select(Variable, Value, DurSym = DelaySym, DurAll = DelayAll))), 
                    cbind(Location = "Kenya", tab_asym_KEN %>% 
                            left_join(tab_dur_KEN %>% select(Variable, Value, DurSym = DelaySym, DurAll = DelayAll)))), 
          file = "docs/tabs/DurAsym.csv")

write.csv(tab_sn_BLT, file = "docs/tabs/DurSn_BLT.csv")
write.csv(tab_sn_KEN, file = "docs/tabs/DurSn_KEN.csv")

write.csv(bind_rows(cbind(Location = "Blantyre", tab_sn_BLT), 
                    cbind(Location = "Kenya", tab_sn_KEN)), 
          file = "docs/tabs/DurSn.csv")


write.csv(tab_sp_BLT, file = "docs/tabs/DurSp_BLT.csv")
write.csv(tab_sp_KEN, file = "docs/tabs/DurSp_KEN.csv")

write.csv(bind_rows(cbind(Location = "Blantyre", tab_sp_BLT), 
                    cbind(Location = "Kenya", tab_sp_KEN)), 
          file = "docs/tabs/DurSp.csv")


write.csv(tab_dur_BLT, file = "docs/tabs/Delay_BLT.csv")
write.csv(tab_dur_KEN, file = "docs/tabs/Delay_KEN.csv")

write.csv(bind_rows(cbind(Location = "Blantyre", tab_dur_BLT), 
                    cbind(Location = "Kenya", tab_dur_KEN)), 
          file = "docs/tabs/Delay.csv")



tab_hiv <- bind_rows(cbind(Location = "Blantyre", tab_dur_BLT), 
                     cbind(Location = "Kenya", tab_dur_KEN)) %>% 
  filter(Variable == "HIV") %>%
  select(Location, Value, PrvAll, CDR_All, IncAll, PrvSym, CDR_Sym, IncSym)
  
t(tab_hiv)[, 4:1]
write.csv(t(tab_hiv)[, 4:1], file = "docs/tabs/HIV.csv")
