
fit_anp <- function(dat_anp, n_iter = 5E4, n_collect = 1E3, n_chain = 3, 
                    prior = list(r_death_sp = 0.127, r_death_sn = 0.024, r_sc_l = 0.1, r_sc_u = 0.3)) {
  require(rstan)
  
  m_anp <- readRDS(file = "stan/m_anp.rds")
  pars <- c("r_sym", "del_sn", "del_sp", "r_tr", "adr", "p_sp", "r_sc")
  
  exo <- prior
  
  ## Duration calc
  f <- "odin/ode_anp_cohort.R"
  model <- odin::odin(f)
  cm <- model()
  calc <- function(cm, inp) {
    cm$set_user(user = inp)
    res <- cm$run(seq(0, 100, 0.1))
    
    a <- sum((res[-1, 2] + res[-nrow(res), 2]) * diff(res[, "t"]) / 2)
    sn <- sum((res[-1, 3] + res[-nrow(res), 3]) * diff(res[, "t"]) / 2)
    sp <- sum((res[-1, 4] + res[-nrow(res), 4]) * diff(res[, "t"]) / 2)
    c(
      A = a,
      Sn = sn,
      Sp = sp,
      Sym = sn + sp,
      All = a + sn + sp, 
      res[nrow(res), c("Cured", "Death", "NotiSn", "NotiSp")]
    )
  }
  
  dat <- c(dat_anp, exo)
  
  fitted <- sampling(m_anp, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain, data = dat)
  mc <- extract(fitted, pars = pars)
  
  sim <- apply(as.data.frame(mc)[c("r_sym", "del_sn", "del_sp", "r_tr", "p_sp", "r_sc")], 1, function(x) {
    calc(cm, c(list(A0 = 1, Sn0 = 0, Sp0 = 0), x))
  })
  sim <- t(sim)
  
  duration <- sim[, c("A", "Sn", "Sp", "Sym", "All")]
  end <- sim[, c("Cured", "Death", "NotiSn", "NotiSp")]
  
  res <- list(
    MC = mc,
    Parameters = summary(fitted, pars = c("r_sym", "r_tr", "adr", "p_sp", "r_sc", "del_sn", "del_sp"))$summary,
    Proportion = summary(fitted, pars = c("pr_a", "pr_sn", "pr_sp"))$summary,
    Prevalence = summary(fitted, pars = "prv")$summary,
    Incidence = list(
      Asym = summary(fitted, pars = "inc_a")$summary,
      Sym = summary(fitted, pars = "inc_s")$summary,
      Sn = summary(fitted, pars = "inc_sn")$summary,
      Sp = summary(fitted, pars = "inc_sp")$summary
    ),
    Notification = list(
      Sn = summary(fitted, pars = "noti_sn")$summary,
      Sp = summary(fitted, pars = "noti_sp")$summary
    ),
    NIR = summary(fitted, pars = "ni")$summary,
    Duration = duration,
    EndPoints = end,
    Data = with(dat, {
      psn <- Sn / N
      psp <- Sp / N
      nsn <- Noti_Sn / Pop
      nsp <- Noti_Sp / Pop
      
      if (year_survey %in% Years) {
        nsn <- nsn[Years == year_survey]
        nsp <- nsp[Years == year_survey]
      } else {
        nsn <- tail(nsn, 1)
        nsp <- tail(nsp, 1)
      }
      
      c(
        PrvSn = psn, PrvSp = psp,
        NotiSn = nsn, NotiSp = nsp,
        DelaySn = psn / nsn, DelaySp = psn / nsn,
        Year = year_survey
      )
    }),
    Exo = prior
  )

  res
}


fit_as <- function(dat_as, n_iter = 5E4, n_collect = 1E3, n_chain = 3, 
                    prior = list(r_death = 0.071, r_sc_l = 0.1, r_sc_u = 0.3)) {
  require(rstan)
  
  m_as <- readRDS(file = "stan/m_as.rds")
  pars <- c("r_sym", "del", "adr", "r_sc")
  
  exo <- prior
  
  ## Duration calc
  f <- "odin/ode_as_cohort.R"
  model <- odin::odin(f)
  cm <- model()
  calc <- function(cm, inp) {
    cm$set_user(user = inp)
    res <- cm$run(seq(0, 100, 0.1))
    
    a <- sum((res[-1, 2] + res[-nrow(res), 2]) * diff(res[, "t"]) / 2)
    sym <- sum((res[-1, 3] + res[-nrow(res), 3]) * diff(res[, "t"]) / 2)

    c(
      A = a,
      Sym = sym,
      All = a + sym, 
      res[nrow(res), c("Cured", "Death", "Noti")]
    )
  }
  
  dat <- c(dat_as, exo)
  
  fitted <- sampling(m_as, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain, data = dat)
  mc <- extract(fitted, pars = c("r_sym", "del", "r_sc", "adr"))
  
  sim <- apply(as.data.frame(mc)[c("r_sym", "del", "r_sc")], 1, function(x) {
    calc(cm, c(list(A0 = 1, S0 = 0), x))
  })
  sim <- t(sim)
  
  duration <- sim[, c("A", "Sym", "All")]
  end <- sim[, c("Cured", "Death", "Noti")]
  
  res <- list(
    MC = mc,
    Parameters = summary(fitted, pars = c("r_sym", "del", "r_sc", "adr"))$summary,
    Proportion = summary(fitted, pars = c("pr_a", "pr_s"))$summary,
    Prevalence = summary(fitted, pars = "prv")$summary,
    Incidence = list(
      Asym = summary(fitted, pars = "inc_a")$summary,
      Sym = summary(fitted, pars = "inc_s")$summary
    ),
    Notification = list(
      Sym = summary(fitted, pars = "noti")$summary
    ),
    #NIR = summary(fitted, pars = "ni")$summary,
    Duration = duration,
    EndPoints = end,
    Data = with(dat, {
      prv <- Sym / N
      noti <- Noti / Pop
      if (year_survey %in% Years) {
        noti <- noti[Years == year_survey]
      } else {
        noti <- tail(noti, 1)
      }
      c(PrvSym = prv, Noti = noti, Delay = prv / noti, Year = year_survey)
    }),
    Exo = prior
  )
  
  res
}


fit_np <- function(dat_np, n_iter = 5E4, n_collect = 1E3, n_chain = 3, 
                   prior = list(r_death_sn = 0.022, r_death_sp = 0.12, r_sc_l = 0.1, r_sc_u = 0.3)) {
  require(rstan)
  
  m_np <- readRDS(file = "stan/m_np.rds")
  pars <- c("del_sn", "del_sp", "r_tr", "adr", "p_sp", "r_sc")
  
  exo <- prior
  
  ## Duration calc
  f <- "odin/ode_np_cohort.R"
  model <- odin::odin(f)
  cm <- model()
  calc <- function(cm, inp) {
    cm$set_user(user = inp)
    res <- cm$run(seq(0, 100, 0.1))
    
    sn <- sum((res[-1, 2] + res[-nrow(res), 2]) * diff(res[, "t"]) / 2)
    sp <- sum((res[-1, 3] + res[-nrow(res), 3]) * diff(res[, "t"]) / 2)
    c(
      Sn = sn,
      Sp = sp,
      Sym = sn + sp,
      All = sn + sp, 
      res[nrow(res), c("Cured", "Death", "NotiSn", "NotiSp")]
    )
  }
  
  dat <- c(dat_np, exo)
  
  fitted <- sampling(m_np, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain, data = dat)
  mc <- extract(fitted, pars = pars)
  
  sim <- apply(as.data.frame(mc)[c("del_sn", "del_sp", "r_tr", "p_sp", "r_sc")], 1, function(x) {
    calc(cm, c(list(n0 = 1), x))
  })
  sim <- t(sim)
  
  duration <- sim[, c("Sn", "Sp", "Sym", "All")]
  end <- sim[, c("Cured", "Death", "NotiSn", "NotiSp")]
  
  res <- list(
    MC = mc,
    Parameters = summary(fitted, pars = c("r_tr", "adr", "p_sp", "r_sc", "del_sn", "del_sp"))$summary,
    Proportion = summary(fitted, pars = c("pr_sn", "pr_sp"))$summary,
    Prevalence = summary(fitted, pars = "prv")$summary,
    Incidence = list(
      Sym = summary(fitted, pars = "inc_s")$summary,
      Sn = summary(fitted, pars = "inc_sn")$summary,
      Sp = summary(fitted, pars = "inc_sp")$summary
    ),
    Notification = list(
      Sn = summary(fitted, pars = "noti_sn")$summary,
      Sp = summary(fitted, pars = "noti_sp")$summary
    ),
    #NIR = summary(fitted, pars = "ni")$summary,
    Duration = duration,
    EndPoints = end,
    Data = with(dat, {
      psn <- Sn / N
      psp <- Sp / N
      nsn <- Noti_Sn / Pop
      nsp <- Noti_Sp / Pop
      
      if (year_survey %in% Years) {
        nsn <- nsn[Years == year_survey]
        nsp <- nsp[Years == year_survey]
      } else {
        nsn <- tail(nsn, 1)
        nsp <- tail(nsp, 1)
      }
      
      c(
        PrvSn = psn, PrvSp = psp,
        NotiSn = nsn, NotiSp = nsp,
        DelaySn = psn / nsn, DelaySp = psn / nsn,
        Year = year_survey
      )
    }),
    Exo = prior
  )
  
  res
}


fit_asc <- function(dat_asc, n_iter = 5E4, n_collect = 1E3, n_chain = 3, 
                   prior = list(r_death = 0.071, r_sc_l = 0.1, r_sc_u = 0.3)) {
  require(rstan)
  
  m_asc <- readRDS(file = "stan/m_asc.rds")
  pars <- c("r_sym", "r_aware", "del", "adr", "r_sc")
  
  exo <- prior
  
  ## Duration calc
  f <- "odin/ode_asc_cohort.R"
  model <- odin::odin(f)
  cm <- model()
  calc <- function(cm, inp) {
    cm$set_user(user = inp)
    res <- cm$run(seq(0, 100, 0.1))
    
    a <- sum((res[-1, 2] + res[-nrow(res), 2]) * diff(res[, "t"]) / 2)
    sym <- sum((res[-1, 3] + res[-nrow(res), 3]) * diff(res[, "t"]) / 2)
    cs <- sum((res[-1, 4] + res[-nrow(res), 4]) * diff(res[, "t"]) / 2)
    
    c(
      A = a,
      Sym = sym,
      CS = cs,
      All = a + sym + cs, 
      res[nrow(res), c("Cured", "Death", "Noti")]
    )
  }
  
  dat <- c(dat_asc, exo)
  
  fitted <- sampling(m_asc, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain, data = dat)
  mc <- extract(fitted, pars = c("r_sym", "r_aware", "del", "r_sc", "adr"))
  
  sim <- apply(as.data.frame(mc)[c("r_sym", "r_aware", "del", "r_sc")], 1, function(x) {
    calc(cm, c(list(A0 = 1, S0 = 0, C0 = 0), x))
  })
  sim <- t(sim)
  
  duration <- sim[, c("A", "Sym", "CS", "All")]
  end <- sim[, c("Cured", "Death", "Noti")]
  
  res <- list(
    MC = mc,
    Parameters = summary(fitted, pars = c("r_sym", "r_aware", "del", "r_sc", "adr"))$summary,
    Proportion = summary(fitted, pars = c("pr_a", "pr_s", "pr_c"))$summary,
    Prevalence = summary(fitted, pars = "prv")$summary,
    Incidence = list(
      Asym = summary(fitted, pars = "inc_a")$summary,
      Sym = summary(fitted, pars = "inc_s")$summary,
      CS = summary(fitted, pars = "inc_c")$summary
    ),
    Notification = list(
      CS = summary(fitted, pars = "noti")$summary
    ),
    #NIR = summary(fitted, pars = "ni")$summary,
    Duration = duration,
    EndPoints = end,
    Data = with(dat, {
      prv <- Sym / N
      noti <- Noti / Pop
      if (year_survey %in% Years) {
        noti <- noti[Years == year_survey]
      } else {
        noti <- tail(noti, 1)
      }
      c(PrvSym = prv, Noti = noti, Delay = prv / noti, Year = year_survey)
    }),
    Exo = prior
  )
  
  res
}



