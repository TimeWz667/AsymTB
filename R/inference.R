source("data/country_list.R")



#i <- 1
#iso <- glue::as_glue(names(countries)[i])
#country <- countries[i]


infer_duration <- function(iso, country) {
  require(tidyverse)
  require(rstan)
  
  ### Load data --
  load("data/Input_" + iso + ".rdata")
  load(paste0("out/ASC/Post_", iso, ".rdata"))
  load(paste0("out/ASC/Post_Sp_", iso, ".rdata"))

  pop <- notification %>%
    filter(Year == 2019) %>%
    group_by(Sex) %>%
    summarise(Pop = sum(Pop)) %>%
    ungroup() %>%
    select(Sex, Pop)
  
  ### Duration ----
  if (country %in% countries_cs) {
    pars <- c("dur_a", "dur_s", "dur_c", "ra", "rs", "rc", "adr", "NotiN", "CDR_A", "CDR_S",
              "IncN_A", "IncN_S", "IncN_C", "PrvN_A", "PrvN_S", "PrvN_C", 
              "ra", "rs", "rc", "Gap_A", "Gap_S", "Gap_C")
    
    pars <- apply(expand.grid(pars, c("[1]", "[2]")), 1, function(x) paste0(x, collapse = ""))
    
    r_sc <- extract(fitted_asc_uni, pars = c("r_sc"))$r_sc
    n_iter <- length(r_sc)
    dur <- as_tibble(extract(fitted_asc_uni, pars = pars)) %>%
      mutate(Key = 1:n_iter, `r_sc[1]` = r_sc, `r_sc[2]` = r_sc) %>%
      pivot_longer(-Key) %>%
      tidyr::extract(name, c("Index", "Sex"), "(\\w+)\\[(\\d)\\]") %>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>%
      pivot_wider(names_from = Index, values_from = value) %>%
      rename(DurA = dur_a, DurS = dur_s, DurC = dur_c) %>%
      mutate(TTS = DurA, TTC = DurA + DurS, TTN = DurA + DurS + DurC,
             Delay_All = TTN, Delay_S = TTN - TTS,
             PrvN_All = PrvN_A + PrvN_S + PrvN_C,
             PrSC_A = DurA * r_sc, PrSC_S = Gap_A * DurS * r_sc, 
             PrSC_C = Gap_A * Gap_S * DurC * r_sc, PrSC = PrSC_A + PrSC_S + PrSC_C,
             PrMor_A = DurA * (ra - r_sc), PrMor_S = Gap_A * DurS * (rs - r_sc),
             PrMor_C = Gap_A * Gap_S * DurC * (rc - r_sc), PrMor = PrMor_A + PrMor_S + PrMor_C
             ) %>%
      select(-c(ra, rs, rc, r_sc, Gap_A, Gap_S, Gap_C))
  } else {
    pars <- c("dur_a", "dur_s", "adr", "NotiN", "CDR_A", "CDR_S", 
              "IncN_A", "IncN_S", "PrvN_A", "PrvN_S",
              "ra", "rs", "Gap_A", "Gap_S")
    
    pars <- apply(expand.grid(pars, c("[1]", "[2]")), 1, function(x) paste0(x, collapse = ""))
    
    r_sc <- extract(fitted_as_uni, pars = c("r_sc"))$r_sc
    n_iter <- length(r_sc)
    dur <- as_tibble(extract(fitted_as_uni, pars = pars)) %>%
      mutate(Key = 1:n_iter, `r_sc[1]` = r_sc, `r_sc[2]` = r_sc) %>%
      pivot_longer(-Key) %>%
      tidyr::extract(name, c("Index", "Sex"), "(\\w+)\\[(\\d)\\]") %>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>%
      pivot_wider(names_from = Index, values_from = value) %>%
      rename(DurA = dur_a, DurS = dur_s) %>%
      mutate(DurC = 0, TTS = DurA, TTC = 0,  TTN = DurA + DurS, 
             Delay_All = TTN, Delay_S = TTN - TTS,
             IncN_C = 0, PrvN_C = 0,
             PrvN_All = PrvN_A + PrvN_S,
             PrSC_A = DurA * r_sc, PrSC_S = Gap_A * DurS * r_sc, 
             PrSC_C = 0, PrSC = PrSC_A + PrSC_S + PrSC_C,
             PrMor_A = DurA * (ra - r_sc), PrMor_S = Gap_A * DurS * (rs - r_sc),
             PrMor_C = 0, PrMor = PrMor_A + PrMor_S + PrMor_C) %>%
      select(-c(ra, rs, r_sc, Gap_A, Gap_S))
  }
  
  dur <- dur %>%
    left_join(pop) %>%
    mutate(Country = country, ISO = iso, 
           Delay_All = TTN, Delay_S = TTN - TTS,
           CNR = NotiN / Pop,
           PN_All = PrvN_All / NotiN, PN_S = (PrvN_S + PrvN_C) / NotiN,
           IncR_All = IncN_A / Pop, IncR_S = IncN_S / Pop,
           IncR_Allg = CNR / CDR_A, IncR_Sg = CNR / CDR_S,
           Prv_All = PrvN_All / Pop) %>%
    rename(ADR = adr)
  
  
  dur_all <- dur %>%
    group_by(Country, ISO, Key) %>%
    summarise(DurA = weighted.mean(DurA , IncN_A), 
              DurS = weighted.mean(DurS, IncN_S), 
              DurC = weighted.mean(DurC, IncN_C),
              ADR = weighted.mean(ADR, IncN_A),
              TTS = weighted.mean(TTS, IncN_A), 
              TTC = weighted.mean(TTC, IncN_A), 
              TTN = weighted.mean(TTN, IncN_A),
              PrMor = weighted.mean(PrMor, IncN_A), 
              PrSC = weighted.mean(PrSC, IncN_A),
              CDR_A = weighted.mean(CDR_A, IncN_A),
              CDR_S = weighted.mean(CDR_S, IncN_S),
              NotiN = sum(NotiN), 
              IncN_A = sum(IncN_A), IncN_S = sum(IncN_S), IncN_C = sum(IncN_C),
              PrvN_A = sum(PrvN_A), PrvN_S = sum(PrvN_S), PrvN_C = sum(PrvN_C), PrvN_All = sum(PrvN_All),
              Pop = sum(Pop)) %>%
    mutate(Delay_All = TTN, Delay_S = TTN - TTS,
           CNR = NotiN / Pop,
           PN_All = PrvN_All / NotiN, PN_S = (PrvN_S + PrvN_C) / NotiN,
           IncR_All = IncN_A / Pop, IncR_S = IncN_S / Pop,
           IncR_Allg = CNR / CDR_A, IncR_Sg = CNR / CDR_S,
           Prv_All = PrvN_All / Pop)
  
  ### Focus on smear-positive -----
  pars <- c("dur_a", "dur_s", "NotiN", "CDR_A", "CDR_S", 
            "IncN_A", "IncN_S", "PrvN_A", "PrvN_S")
  
  pars <- apply(expand.grid(pars, c("[1]", "[2]")), 1, function(x) paste0(x, collapse = ""))
  r_sc <- extract(fitted_as_uni, pars = c("r_sc"))$r_sc
  n_iter <- length(r_sc)
  dur_sp <- as_tibble(extract(fitted_sp, pars = pars)) %>%
    mutate(Key = 1:n_iter) %>%
    pivot_longer(-Key) %>%
    tidyr::extract(name, c("Index", "Sex"), "(\\w+)\\[(\\d)\\]") %>%
    mutate(Sex = ifelse(Sex == 1, "Female", "Male")) %>%
    pivot_wider(names_from = Index, values_from = value) %>%
    left_join(pop) %>%
    rename(Dur_ASp = dur_a, Dur_SymSp = dur_s, PrvN_ASp = PrvN_A, PrvN_SymSp = PrvN_S,
           NotiN_Sp = NotiN, CDR_ASp = CDR_A, CDR_SymSp = CDR_S,
           IncN_ASp = IncN_A, IncN_SymSp = IncN_S) %>%
    mutate(Delay_Sp = Dur_ASp + Dur_SymSp, Delay_SymSp = Dur_SymSp,
           PrvN_Sp = PrvN_ASp + PrvN_SymSp,
           CNR_Sp = NotiN_Sp / Pop,
           PN_Sp = PrvN_Sp / NotiN_Sp, PN_SymSp = PrvN_SymSp / NotiN_Sp)

  dur_sp_all <- dur_sp %>%
    group_by(Key) %>%
    summarise(Pop = sum(Pop),
              Delay_Sp = weighted.mean(Delay_Sp, IncN_ASp), 
              Delay_SymSp = weighted.mean(Delay_SymSp, IncN_SymSp), 
              PrvN_Sp = sum(PrvN_Sp), PrvN_SymSp = sum(PrvN_SymSp), 
              NotiN_Sp = sum(NotiN_Sp), CNR_Sp = sum(NotiN_Sp) / Pop) %>%
    mutate(PN_Sp = PrvN_Sp / NotiN_Sp, PN_SymSp = PrvN_SymSp / NotiN_Sp)
  
  ### Collect results -----
  res <- list()
  
  res$Durations_Sex <- dur %>% left_join(dur_sp)
  
  res$Durations_All <- dur_all %>% left_join(dur_sp_all)
    
  res
}


infer_conversion <- function(iso, country) {
  require(tidyverse)
  require(rstan)
  
  load("out/ANP/Post_" + iso + ".rdata")
  
  res <- list(ISO = iso)
  
  pars <- as_tibble(extract(fitted_anp_uni, c("p_sp", "r_tr")))  
  plot(pars$p_sp, pars$r_tr, main = iso)
  pars$ISO <- iso
  res$Pars <- pars
  
  res$Reg <- lm(p_sp~r_tr, data = pars)
  
  dens <- MASS::kde2d(pars$r_tr, pars$p_sp, n = c(101, 51), lims = c(0, 2, 0, 1))  
  res$Dens <- tibble(x = rep(dens$x, length(dens$y)), 
                     y = rep(dens$y, each = length(dens$x)), 
                     z = c(dens$z),
                     ISO = iso)
  
  res$Wts <- with(dat_anp, {
    c(
      Pop = sum(Pop[, n_t]),
      PrvN = sum(Asym + Sn + Sp) / sum(N) * sum(Pop[, n_t]),
      PrvP = sum(Asym + Sn + Sp) / sum(N),
      NotiN = sum(NotiSn[, n_t] + NotiSp[, n_t]),
      NotiR = sum(NotiSn[, n_t] + NotiSp[, n_t]) / sum(Pop[, n_t])
    )
  })
  
  res
}


infer_cascade <- function(iso, country) {
  require(tidyverse)
  require(rstan)
  
  ### Load data -----
  load("out/ASC/Post_" + iso + ".rdata")
  load("data/TO_" + iso + ".rdata")
  
  gap_to <- tibble(Sex = c("Female", "Male"),
                   Gap_T = pars_to$r_cure / sum(unlist(pars_to)))
  
  if (country %in% countries_cs) {
    wts <- extract(fitted_asc_uni, pars = c("IncN_A"))$IncN_A
    n_iter <- dim(wts)[1]
    
    wts <- tibble(Key = rep(1:n_iter, 2), 
                  Sex = rep(c("Female", "Male"), each = n_iter), 
                  Wts = c(wts / rowSums(wts)))
    
    pars <- c("Gap_A", "Gap_S", "Gap_C")
    pars <- apply(expand.grid(pars, c("[1]", "[2]")), 1, function(x) paste0(x, collapse = ""))
    
    
    Cascade_Sex <- as_tibble(extract(fitted_asc_uni, pars = pars)) %>%
      mutate(Key = 1:n_iter) %>%
      pivot_longer(-Key) %>%
      tidyr::extract(name, c("Index", "Sex"), "(\\w+)\\[(\\d)\\]") %>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male"), Country = country, ISO = iso) %>%
      pivot_wider(names_from = Index, values_from = value) %>%
      left_join(gap_to) %>%
      rename(G_sym = Gap_A, G_aware = Gap_S, G_care = Gap_C, G_complete = Gap_T) %>%
      mutate(C_inc = 1, G_inc = 1, 
             C_sym = G_sym, C_aware = G_aware * C_sym, 
             C_care = G_care * C_aware, C_complete = G_complete * C_care,) %>%
      left_join(wts)
  } else {
    wts <- extract(fitted_as_uni, pars = c("IncN_A"))$IncN_A
    n_iter <- dim(wts)[1]
    
    wts <- tibble(Key = rep(1:n_iter, 2), 
                  Sex = rep(c("Female", "Male"), each = n_iter), 
                  Wts = c(wts / rowSums(wts)))
    
    pars <- c("Gap_A", "Gap_S")
    pars <- apply(expand.grid(pars, c("[1]", "[2]")), 1, function(x) paste0(x, collapse = ""))
    
    Cascade_Sex <- as_tibble(extract(fitted_as_uni, pars = pars)) %>%
      mutate(Key = 1:n_iter) %>%
      pivot_longer(-Key) %>%
      tidyr::extract(name, c("Index", "Sex"), "(\\w+)\\[(\\d)\\]") %>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male"), Country = country, ISO = iso) %>%
      pivot_wider(names_from = Index, values_from = value) %>%
      left_join(gap_to) %>%
      rename(G_sym = Gap_A, G_care = Gap_S, G_complete = Gap_T) %>%
      mutate(C_inc = 1, G_inc = 1, 
             C_sym = G_sym, C_care = G_care * C_sym, C_complete = G_complete * C_care,
             G_aware = 0, C_aware = 0) %>%
      left_join(wts)
  }


  Cascade_All <- Cascade_Sex %>%
    group_by(Key, Country, ISO) %>%
    summarise(
      G_inc = sum(G_inc * Wts), 
      G_sym = sum(G_sym * Wts), G_aware = sum(G_aware * Wts), 
      G_care = sum(G_care * Wts), G_complete = sum(G_complete * Wts),
      C_inc = sum(C_inc * Wts), 
      C_sym = sum(C_sym * Wts), C_aware = sum(C_aware * Wts), 
      C_care = sum(C_care * Wts), C_complete = sum(C_complete * Wts)
    )
  
  list(
    Cascade_All = Cascade_All,
    Cascade_Sex = Cascade_Sex
  )
}


infer_cohort <- function(iso, country, n_sim = 200) {
  require(odin)
  iso <- glue::as_glue(iso)

  
  ### Load data ----
  load("out/ASC/Post_" + iso + ".rdata")
  load("data/TO_" + iso + ".rdata")
  
  if (country %in% countries_cs) {
    wts <- extract(fitted_asc_uni, pars = "IncN_A")$IncN_A
    wts <- wts / rowSums(wts)
  } else {
    wts <- extract(fitted_as_uni, pars = "IncN_A")$IncN_A
    wts <- wts / rowSums(wts)
  }
  
  
  
  ### Cascade, time ----
  if (country %in% countries_cs) {
    
    pars <- extract(fitted_asc_uni, c("r_sym", "r_aware", "r_det", "r_sc"))

    f <- "odin/ode_asct_cohort.R"
    model <- odin(f)
    
    ys_ts <- data.table::rbindlist(lapply(1:2, function(sex) {
      sims <- lapply(1:n_sim, function(j) {
        inp <- list(
          r_sym = pars$r_sym[j, sex],
          r_aware = pars$r_aware[j, sex],
          r_sc = pars$r_sc[j, sex],
          r_det = pars$r_det[j, sex],
          r_cure = pars_to$r_cure,
          r_ltfu = pars_to$r_ltfu,
          r_death_a = dat_asc$r_death_a[sex],
          r_death_s = dat_asc$r_death_s[sex],
          r_death_ontr = pars_to$r_tbmu_ontr
        )
        cm <- model(user = inp)
        
        ys <- cm$run(seq(0, 2, 0.02))
        
        data.table::data.table(ys) %>%
          pivot_longer(cols = A:Cured)
      })
      
      data.table::rbindlist(sims) %>%
        group_by(t, name) %>%
        summarise(value = mean(value)) %>%
        mutate(Sex = sex)
      
    })) %>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male"),
             name = factor(name, levels = c("Death", "LTFU", "SelfCured", "Cured", "Tr", "C", "S", "A")))
    
    ys_ts_total <- ys_ts %>% left_join(data.frame(Sex = c("Female", "Male"), wts = colMeans(wts))) %>%
      group_by(t, name) %>%
      summarise(value = weighted.mean(value, wts)) %>%
      mutate(Sex = "Total")
    
    
    ys_end <- data.table::rbindlist(lapply(1:2, function(sex) {
      sims <- lapply(1:n_sim, function(j) {
        inp <- list(
          r_sym = pars$r_sym[j, sex],
          r_aware = pars$r_aware[j, sex],
          r_sc = pars$r_sc[j, sex],
          r_det = pars$r_det[j, sex],
          r_cure = pars_to$r_cure,
          r_ltfu = pars_to$r_ltfu,
          r_death_a = dat_asc$r_death_a[sex],
          r_death_s = dat_asc$r_death_s[sex],
          r_death_ontr = pars_to$r_tbmu_ontr
        )
        cm <- model(user = inp)
        
        ys <- cm$run(seq(0, 20, 0.5))
        
        data.table::data.table(ys) %>%
          pivot_longer(cols = A:Cured) %>%
          filter(t == 20)
      })
      
      data.table::rbindlist(sims) %>%
        group_by(t, name) %>%
        summarise(value = mean(value)) %>%
        mutate(Sex = sex)
      
    })) %>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male"),
             name = factor(name, levels = c("Death", "LTFU", "SelfCured", "Cured", "Tr", "C", "S", "A"))) %>%
      filter(name %in% c("Death", "LTFU", "SelfCured", "Cured"))
    
    
    ys_end_total <- ys_end %>% left_join(data.frame(Sex = c("Female", "Male"), wts = colMeans(wts))) %>%
      group_by(t, name) %>%
      summarise(value = weighted.mean(value, wts)) %>%
      mutate(Sex = "Total")
  } else {
    pars <- extract(fitted_as_uni, c("r_sym", "r_det", "r_sc"))
    f <- "odin/ode_ast_cohort.R"
    model <- odin(f)
    
    ys_ts <- data.table::rbindlist(lapply(1:2, function(sex) {
      sims <- lapply(1:n_sim, function(j) {
        inp <- list(
          r_sym = pars$r_sym[j, sex],
          r_sc = pars$r_sc[j, sex],
          r_det = pars$r_det[j, sex],
          r_cure = pars_to$r_cure,
          r_ltfu = pars_to$r_ltfu,
          r_death_a = dat_as$r_death_a[sex],
          r_death_s = dat_as$r_death_s[sex],
          r_death_ontr = pars_to$r_tbmu_ontr
        )
        cm <- model(user = inp)
        
        ys <- cm$run(seq(0, 2, 0.02))
        
        data.table::data.table(ys) %>%
          pivot_longer(cols = A:Cured)
      })
      
      data.table::rbindlist(sims) %>%
        group_by(t, name) %>%
        summarise(value = mean(value)) %>%
        mutate(Sex = sex)
      
    })) %>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male"),
             name = factor(name, levels = c("Death", "LTFU", "SelfCured", "Cured", "Tr", "S", "A")))
    
    ys_ts_total <- ys_ts %>% left_join(data.frame(Sex = c("Female", "Male"), wts = colMeans(wts))) %>%
      group_by(t, name) %>%
      summarise(value = weighted.mean(value, wts)) %>%
      mutate(Sex = "Total")
    
    
    ys_end <- data.table::rbindlist(lapply(1:2, function(sex) {
      sims <- lapply(1:n_sim, function(j) {
        inp <- list(
          r_sym = pars$r_sym[j, sex],
          r_sc = pars$r_sc[j, sex],
          r_det = pars$r_det[j, sex],
          r_cure = pars_to$r_cure,
          r_ltfu = pars_to$r_ltfu,
          r_death_a = dat_as$r_death_a[sex],
          r_death_s = dat_as$r_death_s[sex],
          r_death_ontr = pars_to$r_tbmu_ontr
        )
        cm <- model(user = inp)
        
        ys <- cm$run(seq(0, 20, 0.5))
        
        data.table::data.table(ys) %>%
          pivot_longer(cols = A:Cured) %>%
          filter(t == 20)
      })
      
      data.table::rbindlist(sims) %>%
        group_by(t, name) %>%
        summarise(value = mean(value)) %>%
        mutate(Sex = sex)
      
    })) %>%
      mutate(Sex = ifelse(Sex == 1, "Female", "Male"),
             name = factor(name, levels = c("Death", "LTFU", "SelfCured", "Cured", "Tr", "S", "A"))) %>%
      filter(name %in% c("Death", "LTFU", "SelfCured", "Cured"))
    
    
    ys_end_total <- ys_end %>% left_join(data.frame(Sex = c("Female", "Male"), wts = colMeans(wts))) %>%
      group_by(t, name) %>%
      summarise(value = weighted.mean(value, wts)) %>%
      mutate(Sex = "Total")
  }
  
  list(
    Cohort = rbind(ys_ts, ys_ts_total) %>% mutate(Country = country, ISO = iso),
    End = rbind(ys_end, ys_end_total) %>% mutate(Country = country, ISO = iso)
  )
}
