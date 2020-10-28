library(rstan)




sub <- glue::as_glue("KEN")


fetch_dur_a <- function(sub) {
  frm_mci <- function(ds, digits = 1) {
    sprintf("%s (%s-%s)", 
            format(mean(ds, na.rm = T), nsmall = digits, digits = digits), 
            format(quantile(ds, 0.025, na.rm = T), nsmall = digits, digits = digits), 
            format(quantile(ds, 0.975, na.rm = T), nsmall = digits, digits = digits))
  }
  
  
  
  sub <- glue::as_glue(sub)
  
  load(file = "out/Full/" + sub +"/Cov.rdata")

  load(file = "out/Full/" + sub +"/Sex.rdata")
  
  ext <- extract(fitted1, pars = c("r_sym", "dur_a"))
  exo <- dataset$exo
  

  tab_sex <- data.table::data.table(
    Variable = "Sex",
    Value = dataset$prv$Sex, 
    Duration = apply(ext$dur_a * 12, 2, frm_mci),
    Rate = apply(ext$r_sym, 2, frm_mci, digit = 2),
    RR_uni = c("reference", frm_mci(exp(extract(fitted3, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 1]), digit = 2))
  )
  tab_sex
  
  
  load(file = "out/Full/" + sub +"/Age.rdata")
  
  ext <- extract(fitted1, pars = c("r_sym", "dur_a"))
  exo <- dataset$exo
  
  n_agp <- length(unique(dataset$prv$Agp))
  if (n_agp > 2) {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = dataset$prv$Agp, 
      Duration = apply(ext$dur_a * 12, 2, frm_mci),
      Rate = apply(ext$r_sym, 2, frm_mci, digit = 2),
      RR_uni = c("reference", apply(exp(extract(fitted3, pars = c("lrr_sym"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2)),
      RR_multi = c("reference", apply(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2))
    )
  } else {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = dataset$prv$Agp, 
      Duration = apply(ext$dur_a * 12, 2, frm_mci),
      Rate = apply(ext$r_sym, 2, frm_mci, digit = 2),
      RR_uni = c("reference", frm_mci(exp(extract(fitted3, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
      RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 2]), digit = 2))
    )
  }

  
  load(file = "out/Full/" + sub + "/HIV.rdata")
  
  ext <- extract(fitted1, pars = c("r_sym", "dur_a"))
  exo <- dataset$exo
  
  tab_hiv <- data.table::data.table(
    Variable = "HIV",
    Value = rev(dataset$prv$HIV), 
    Duration = rev(apply(ext$dur_a * 12, 2, frm_mci)),
    Rate = rev(apply(ext$r_sym, 2, frm_mci, digit = 2)),
    RR_uni = c("reference", frm_mci(1/exp(extract(fitted3, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, n_agp + 1]), digit = 2))
  )
  tab_hiv
  
  
  load(file = "out/Full/" + sub + "/Total.rdata")
  ext <- extract(fitted3, pars = c("r_sym", "dur_a"))
  
  
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
  frm_mci <- function(ds, digits = 1) {
    sprintf("%s (%s-%s)", 
            format(mean(ds, na.rm = T), nsmall = digits, digits = digits), 
            format(quantile(ds, 0.025, na.rm = T), nsmall = digits, digits = digits), 
            format(quantile(ds, 0.975, na.rm = T), nsmall = digits, digits = digits))
  }
  
  
  sub <- glue::as_glue(sub)
  
  load(file = "out/Full/" + sub +"/Cov.rdata")
  
  load(file = "out/Full/" + sub +"/Sex.rdata")
  
  ext <- extract(fitted1, pars = c("r_det_sn", "dur_sn"))
  exo <- dataset$exo
  
  
  tab_sex <- data.table::data.table(
    Variable = "Sex",
    Value = dataset$prv$Sex, 
    Duration = apply(ext$dur_sn * 12, 2, frm_mci),
    Rate = apply(ext$r_det_sn, 2, frm_mci, digit = 2),
    RR_uni = c("reference", frm_mci(exp(extract(fitted3, pars = c("lrr_cs_sn"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sn"))[[1]][, 1]), digit = 2))
  )
  tab_sex
  
  
  load(file = "out/Full/" + sub +"/Age.rdata")
  
  ext <- extract(fitted1, pars = c("r_det_sn", "dur_sn"))
  exo <- dataset$exo
  
  n_agp <- length(unique(dataset$prv$Agp))
  if (n_agp > 2) {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = dataset$prv$Agp, 
      Duration = apply(ext$dur_sn * 12, 2, frm_mci),
      Rate = apply(ext$r_det_sn, 2, frm_mci, digit = 2),
      RR_uni = c("reference", apply(exp(extract(fitted3, pars = c("lrr_cs_sn"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2)),
      RR_multi = c("reference", apply(exp(extract(fitted_full, c("lrr_cs_sn"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2))
    )
  } else {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = dataset$prv$Agp, 
      Duration = apply(ext$dur_sn * 12, 2, frm_mci),
      Rate = apply(ext$r_det_sn, 2, frm_mci, digit = 2),
      RR_uni = c("reference", frm_mci(exp(extract(fitted3, pars = c("lrr_cs_sn"))[[1]][, 2]), digit = 2)),
      RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sn"))[[1]][, 2]), digit = 2))
    )
  }
  
  
  load(file = "out/Full/" + sub + "/HIV.rdata")
  
  ext <- extract(fitted1, pars = c("r_det_sn", "dur_sn"))
  exo <- dataset$exo
  
  tab_hiv <- data.table::data.table(
    Variable = "HIV",
    Value = rev(dataset$prv$HIV), 
    Duration = rev(apply(ext$dur_sn * 12, 2, frm_mci)),
    Rate = rev(apply(ext$r_det_sn, 2, frm_mci, digit = 2)),
    RR_uni = c("reference", frm_mci(1/exp(extract(fitted3, pars = c("lrr_cs_sn"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sn"))[[1]][, n_agp + 1]), digit = 2))
  )
  tab_hiv
  
  
  load(file = "out/Full/" + sub + "/Total.rdata")
  ext <- extract(fitted3, pars = c("r_det_sn", "dur_sn"))
  
  
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
  frm_mci <- function(ds, digits = 1) {
    sprintf("%s (%s-%s)", 
            format(mean(ds, na.rm = T), nsmall = digits, digits = digits), 
            format(quantile(ds, 0.025, na.rm = T), nsmall = digits, digits = digits), 
            format(quantile(ds, 0.975, na.rm = T), nsmall = digits, digits = digits))
  }
  
  
  sub <- glue::as_glue(sub)
  
  load(file = "out/Full/" + sub +"/Cov.rdata")
  
  load(file = "out/Full/" + sub +"/Sex.rdata")
  
  ext <- extract(fitted1, pars = c("r_det_sp", "dur_sp"))
  exo <- dataset$exo
  
  
  tab_sex <- data.table::data.table(
    Variable = "Sex",
    Value = dataset$prv$Sex, 
    Duration = apply(ext$dur_sp * 12, 2, frm_mci),
    Rate = apply(ext$r_det_sp, 2, frm_mci, digit = 2),
    RR_uni = c("reference", frm_mci(exp(extract(fitted3, pars = c("lrr_cs_sp"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sp"))[[1]][, 1]), digit = 2))
  )
  tab_sex
  
  
  load(file = "out/Full/" + sub +"/Age.rdata")
  
  ext <- extract(fitted1, pars = c("r_det_sp", "dur_sp"))
  exo <- dataset$exo
  
  n_agp <- length(unique(dataset$prv$Agp))
  if (n_agp > 2) {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = dataset$prv$Agp, 
      Duration = apply(ext$dur_sp * 12, 2, frm_mci),
      Rate = apply(ext$r_det_sp, 2, frm_mci, digit = 2),
      RR_uni = c("reference", apply(exp(extract(fitted3, pars = c("lrr_cs_sp"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2)),
      RR_multi = c("reference", apply(exp(extract(fitted_full, c("lrr_cs_sp"))[[1]][, 2:n_agp]), 2, frm_mci, digit = 2))
    )
  } else {
    tab_age <- data.table::data.table(
      Variable = "Age",
      Value = dataset$prv$Agp, 
      Duration = apply(ext$dur_sp * 12, 2, frm_mci),
      Rate = apply(ext$r_det_sp, 2, frm_mci, digit = 2),
      RR_uni = c("reference", frm_mci(exp(extract(fitted3, pars = c("lrr_cs_sp"))[[1]][, 2]), digit = 2)),
      RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sp"))[[1]][, 2]), digit = 2))
    )
  }
  
  
  load(file = "out/Full/" + sub + "/HIV.rdata")
  
  ext <- extract(fitted1, pars = c("r_det_sp", "dur_sp"))
  exo <- dataset$exo
  
  tab_hiv <- data.table::data.table(
    Variable = "HIV",
    Value = rev(dataset$prv$HIV), 
    Duration = rev(apply(ext$dur_sp * 12, 2, frm_mci)),
    Rate = rev(apply(ext$r_det_sp, 2, frm_mci, digit = 2)),
    RR_uni = c("reference", frm_mci(1/exp(extract(fitted3, pars = c("lrr_cs_sp"))[[1]][, 2]), digit = 2)),
    RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_cs_sp"))[[1]][, n_agp + 1]), digit = 2))
  )
  tab_hiv
  
  
  load(file = "out/Full/" + sub + "/Total.rdata")
  ext <- extract(fitted3, pars = c("r_det_sp", "dur_sp"))
  
  
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



tab_asym_KEN <- fetch_dur_a("KEN")
tab_sn_KEN <- fetch_dur_sn("KEN")
tab_sp_KEN <- fetch_dur_sp("KEN")


tab_asym_BLT <- fetch_dur_a("BLT")
tab_sn_BLT <- fetch_dur_sn("BLT")
tab_sp_BLT <- fetch_dur_sp("BLT")




#### Output ----
write.csv(tab_asym_BLT, file = "docs/tabs/DurAsym_BLT.csv")
write.csv(tab_asym_KEN, file = "docs/tabs/DurAsym_KEN.csv")

write.csv(rbind(cbind(Location = "Blantyre", tab_asym_BLT), 
                cbind(Location = "Kenya", tab_asym_KEN)), 
          file = "docs/tabs/DurAsym.csv")

write.csv(tab_sn_BLT, file = "docs/tabs/DurSn_BLT.csv")
write.csv(tab_sn_KEN, file = "docs/tabs/DurSn_KEN.csv")

write.csv(rbind(cbind(Location = "Blantyre", tab_sn_BLT), 
                cbind(Location = "Kenya", tab_sn_KEN)), 
          file = "docs/tabs/DurSn.csv")


write.csv(tab_sp_BLT, file = "docs/tabs/DurSp_BLT.csv")
write.csv(tab_sp_KEN, file = "docs/tabs/DurSp_KEN.csv")

write.csv(rbind(cbind(Location = "Blantyre", tab_sp_BLT), 
                cbind(Location = "Kenya", tab_sp_KEN)), 
          file = "docs/tabs/DurSp.csv")

