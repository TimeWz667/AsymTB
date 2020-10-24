library(rstan)


frm_mci <- function(ds, digits = 1) {
  sprintf("%s (%s-%s)", 
          format(mean(ds, na.rm = T), nsmall = digits, digits = digits), 
          format(quantile(ds, 0.025, na.rm = T), nsmall = digits, digits = digits), 
          format(quantile(ds, 0.975, na.rm = T), nsmall = digits, digits = digits))
}

#### Kenya ----
load(file = "out/Full/KEN/Cov.rdata")


load(file = "out/Full/KEN/Sex.rdata")

ext <- extract(fitted1, pars = c("r_sym", "r_sc"))
exo <- dataset$exo

durations <- sapply(1:nrow(ext$r_sym), function(i) {
  r_sym <- ext$r_sym[i, ]
  r_sc <- ext$r_sc[i]
  r_death_bg <- exo$r_death_bg
  12 / (r_sym + r_sc + r_death_bg)
})

tab_sex <- data.table::data.table(
  Variable = "Sex",
  Value = dataset$prv$Sex, 
  Duration = apply(durations, 1, frm_mci),
  Rate = apply(ext$r_sym, 2, frm_mci, digit = 2),
  RR_uni = c("reference", frm_mci(exp(extract(fitted3, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
  RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 1]), digit = 2))
)
tab_sex



load(file = "out/Full/KEN/Age.rdata")

ext <- extract(fitted1, pars = c("r_sym", "r_sc"))
exo <- dataset$exo

durations <- sapply(1:nrow(ext$r_sym), function(i) {
  r_sym <- ext$r_sym[i, ]
  r_sc <- ext$r_sc[i]
  r_death_bg <- exo$r_death_bg
  12 / (r_sym + r_sc + r_death_bg)
})

tab_age <- data.table::data.table(
  Variable = "Age",
  Value = dataset$prv$Agp, 
  Duration = apply(durations, 1, frm_mci),
  Rate = apply(ext$r_sym, 2, frm_mci, digit = 2),
  RR_uni = c("reference", apply(exp(extract(fitted3, pars = c("lrr_sym"))[[1]][, 2:5]), 2, frm_mci, digit = 2)),
  RR_multi = c("reference", apply(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 2:5]), 2, frm_mci, digit = 2))
)
tab_age


load(file = "out/Full/KEN/HIV.rdata")

ext <- extract(fitted1, pars = c("r_sym", "r_sc"))
exo <- dataset$exo

durations <- sapply(1:nrow(ext$r_sym), function(i) {
  r_sym <- ext$r_sym[i, ]
  r_sc <- ext$r_sc[i]
  r_death_bg <- exo$r_death_bg
  12 / (r_sym + r_sc + r_death_bg)
})

tab_hiv <- data.table::data.table(
  Variable = "HIV",
  Value = rev(dataset$prv$HIV), 
  Duration = rev(apply(durations, 1, frm_mci)),
  Rate = rev(apply(ext$r_sym, 2, frm_mci, digit = 2)),
  RR_uni = c("reference", frm_mci(1/exp(extract(fitted3, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
  RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 6]), digit = 2))
)
tab_hiv


load(file = "output/Full/KEN/Total.rdata")
ext <- extract(fitted3, pars = c("r_sym", "r_sc"))

durations <- sapply(1:length(ext$r_sym), function(i) {
  r_sym <- ext$r_sym[i]
  r_sc <- ext$r_sc[i]
  r_death_bg <- exo$r_death_bg
  12 / (r_sym + r_sc + r_death_bg)
})


tab_overall <- data.table::data.table(
  Variable = "Overall",
  Value = "", 
  Duration = frm_mci(durations),
  Rate = frm_mci(ext$r_sym, digit = 2),
  RR_uni = "",
  RR_multi = ""
)

tab_kenya <- rbind(tab_overall, tab_sex, tab_age, tab_hiv)



#### Blantyre ----
load(file = "out/Full/BLT/Cov.rdata")


load(file = "out/Full/BLT/Sex.rdata")

ext <- extract(fitted1, pars = c("r_sym", "r_sc"))
exo <- dataset$exo

durations <- sapply(1:nrow(ext$r_sym), function(i) {
  r_sym <- ext$r_sym[i, ]
  r_sc <- ext$r_sc[i]
  r_death_bg <- exo$r_death_bg
  12 / (r_sym + r_sc + r_death_bg)
})

tab_sex <- data.table::data.table(
  Variable = "Sex",
  Value = dataset$prv$Sex, 
  Duration = apply(durations, 1, frm_mci),
  Rate = apply(ext$r_sym, 2, frm_mci, digit = 2),
  RR_uni = c("reference", frm_mci(exp(extract(fitted3, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
  RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 1]), digit = 2))
)
tab_sex



load(file = "out/Full/BLT/Age.rdata")

ext <- extract(fitted1, pars = c("r_sym", "r_sc"))
exo <- dataset$exo

durations <- sapply(1:nrow(ext$r_sym), function(i) {
  r_sym <- ext$r_sym[i, ]
  r_sc <- ext$r_sc[i]
  r_death_bg <- exo$r_death_bg
  12 / (r_sym + r_sc + r_death_bg)
})

tab_age <- data.table::data.table(
  Variable = "Age",
  Value = dataset$prv$Agp, 
  Duration = apply(durations, 1, frm_mci),
  Rate = apply(ext$r_sym, 2, frm_mci, digit = 2),
  RR_uni = c("reference", frm_mci(exp(extract(fitted3, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
  RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 2]), digit = 2))
)
tab_age


load(file = "out/Full/BLT/HIV.rdata")

ext <- extract(fitted1, pars = c("r_sym", "r_sc"))
exo <- dataset$exo

durations <- sapply(1:nrow(ext$r_sym), function(i) {
  r_sym <- ext$r_sym[i, ]
  r_sc <- ext$r_sc[i]
  r_death_bg <- exo$r_death_bg
  12 / (r_sym + r_sc + r_death_bg)
})

tab_hiv <- data.table::data.table(
  Variable = "HIV",
  Value = rev(dataset$prv$HIV), 
  Duration = rev(apply(durations, 1, frm_mci)),
  Rate = rev(apply(ext$r_sym, 2, frm_mci, digit = 2)),
  RR_uni = c("reference", frm_mci(1/exp(extract(fitted3, pars = c("lrr_sym"))[[1]][, 2]), digit = 2)),
  RR_multi = c("reference", frm_mci(exp(extract(fitted_full, c("lrr_sym"))[[1]][, 3]), digit = 2))
)
tab_hiv


load(file = "out/Full/BLT/Total.rdata")
ext <- extract(fitted3, pars = c("r_sym", "r_sc"))

durations <- sapply(1:length(ext$r_sym), function(i) {
  r_sym <- ext$r_sym[i]
  r_sc <- ext$r_sc[i]
  r_death_bg <- exo$r_death_bg
  12 / (r_sym + r_sc + r_death_bg)
})


tab_overall <- data.table::data.table(
  Variable = "Overall",
  Value = "", 
  Duration = frm_mci(durations),
  Rate = frm_mci(ext$r_sym, digit = 2),
  RR_uni = "",
  RR_multi = ""
)

tab_blantyre <- rbind(tab_overall, tab_sex, tab_age, tab_hiv)



#### Output ----
write.csv(tab_blantyre, file = "docs/tabs/DurAsym_BLT.csv")
write.csv(tab_kenya, file = "docs/tabs/DurAsym_KEN.csv")

write.csv(rbind(cbind(Location = "Blantyre", tab_blantyre), 
                cbind(Location = "Kenya", tab_kenya)), 
          file = "docs/tabs/DurAsym.csv")
