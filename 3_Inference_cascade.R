library(odin)
library(tidyverse)
library(rstan)

### Settings -----
options(odin.verbose = F,
        odin.target = ifelse(odin::can_compile(), "c", "r"),
        odin.compiler_warnings = F,
        odin.no_check_unused_equations = T
)

f <- "odin/ode_ast_cohort.R"
m_ast <- odin(f)

f <- "odin/ode_asct_cohort.R"
m_asct <- odin(f)

n_sim <- 100


### Country list -----
countries <- c(
  KHM = "Cambodia",
  KEN = "Kenya",
  LAO = "Lao People's Democratic Republic", 
  MWI = "Malawi", 
  PAK = "Pakistan", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda", #
  VNM = "Viet Nam", 
  ZMB = "Zambia"
)

countries_cs <- c(
  KEN = "Kenya",
  MWI = "Malawi", 
  PHL = "Philippines", 
  TZA = "United Republic of Tanzania", 
  UGA = "Uganda", 
  ZMB = "Zambia"
)





for (i in 1:length(countries)) {
  iso <- glue::as_glue(names(countries)[i])
  country <- countries[i]
  
  ### Load data ----
  load("out/ASC/Post_" + iso + ".rdata")
  load("data/TO_" + iso + ".rdata")
  
  if (country %in% countries_cs) {
    wts <- extract(fitted_asc_uni, pars = "inc_a")$inc_a
    wts <- wts[, , dim(wts)[3]]
    wts <- wts / rowSums(wts)
  } else {
    wts <- extract(fitted_as_uni, pars = "inc_a")$inc_a
    wts <- wts[, , dim(wts)[3]]
    wts <- wts / rowSums(wts)
  }
  
  
  
  ### Cascade, time ----
  if (country %in% countries_cs) {
 
    pars <- extract(fitted_asc_uni, c("r_sym", "r_aware", "r_det", "r_sc"))
    
    model <- m_asct
    
    ys_ts <- data.table::rbindlist(lapply(1:2, function(sex) {
      sims <- lapply(1:n_sim, function(j) {
        inp <- list(
          r_sym = pars$r_sym[j, sex],
          r_aware = pars$r_aware[j, sex],
          r_sc = pars$r_sc[j, sex],
          r_det = pars$r_det[j, sex],
          r_cure = pars_to$r_cure,
          r_ltfu = pars_to$r_ltfu,
          r_death_bg = dat_asc$r_death_bg[sex],
          r_death_untr = dat_asc$r_death_tb[sex],
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
          r_death_bg = dat_asc$r_death_bg[sex],
          r_death_untr = dat_asc$r_death_tb[sex],
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
    
    model <- m_ast
    
    ys_ts <- data.table::rbindlist(lapply(1:2, function(sex) {
      sims <- lapply(1:n_sim, function(j) {
        inp <- list(
          r_sym = pars$r_sym[j, sex],
          r_sc = pars$r_sc[j, sex],
          r_det = pars$r_det[j, sex],
          r_cure = pars_to$r_cure,
          r_ltfu = pars_to$r_ltfu,
          r_death_bg = dat_as$r_death_bg[sex],
          r_death_untr = dat_as$r_death_tb[sex],
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
          r_death_bg = dat_as$r_death_bg[sex],
          r_death_untr = dat_as$r_death_tb[sex],
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
  
  Cohort <- list(
    TS = rbind(ys_ts, ys_ts_total),
    End = rbind(ys_end, ys_end_total)
  )
  ### Cascade, percent ----
  
  if (country %in% countries_cs) {
    ext <- extract(fitted_asc_uni, pars = c("ra", "r_sym", "rs", "r_aware", "rc", "r_det"))
    
    gaps <- with(ext, {
      data.frame(G_sym = rowSums(r_sym / (ra + r_sym) * wts), 
                 G_aware = rowSums(r_aware / (rs + r_aware) * wts), 
                 G_care = rowSums(r_det / (rc + r_det) * wts),
                 G_complete = pars_to$r_cure / sum(unlist(pars_to))
      )
    })
    
    cas <- t(apply(gaps, 1, cumprod))
    colnames(cas) <- c("C_sym", "C_aware", "C_care", "C_complete")
    
    
    Cascade <- cbind(gaps, cas) %>%
      summarise(across(everything(), list(M = mean, 
                                          L = ~quantile(., 0.025),
                                          U = ~quantile(., 0.975)))) %>%
      pivot_longer(starts_with(c("C_", "G_")), values_to = "pr") %>%
      separate(name, c("Index", "Stage", "Stat"), "_") %>%
      pivot_wider(c(Index, Stage), names_from = Stat, values_from = pr) %>%
      mutate(Stage = factor(Stage, levels = c("sym", "aware", "care", "complete")))
    
    
  } else {
    ext <- extract(fitted_as_uni, pars = c("ra", "r_sym", "rs", "r_det"))
    
    gaps <- with(ext, {
      data.frame(G_sym = rowSums(r_sym / (ra + r_sym) * wts), 
                 G_care = rowSums(r_det / (rs + r_det) * wts),
                 G_complete = pars_to$r_cure / sum(unlist(pars_to))
      )
    })
    
    cas <- t(apply(gaps, 1, cumprod))
    colnames(cas) <- c("C_sym", "C_care", "C_complete")
    
    
    Cascade <- cbind(gaps, cas) %>%
      summarise(across(everything(), list(M = mean, 
                                          L = ~quantile(., 0.025),
                                          U = ~quantile(., 0.975)))) %>%
      pivot_longer(starts_with(c("C_", "G_")), values_to = "pr") %>%
      separate(name, c("Index", "Stage", "Stat"), "_") %>%
      pivot_wider(c(Index, Stage), names_from = Stat, values_from = pr) %>%
      mutate(Stage = factor(Stage, levels = c("sym", "aware", "care", "complete")))
  }
  
  save(Cascade, Cohort, file = "out/Cascade_" + iso + ".rdata")
}





