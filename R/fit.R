## Compile models -----
# m_at_ts <- stan_model("stan/m_anp_ts.stan")
# saveRDS(m_at_ts, file = "stan/m_at_ts.rds")
# m_t_ts <- stan_model("stan/m_np_ts.stan")
# saveRDS(m_t_ts, file = "stan/m_t_ts.rds")


## Model fitting
fit_anp <- function(dat_anp, dur_asym = 3, n_iter = 5E4, n_collect = 1E3, n_chain = 3) {
  require(rstan)
  
  m_at_ts <- readRDS(file = "stan/m_at_ts.rds")
  
  
  pars <- c("del_sp", "del_sn", "r_tr", "adr", "p_sp", "pr_sn", "pr_sp", "pr_a",
            "inc_a", "inc_s", "inc_sn", "inc_sp", "noti_sn", "noti_sp", "prv")
  
  
  inp <- list(r_death_sp = 0.127, r_death_sn = 0.024, r_sym = 12/dur_asym, r_sc = 0.15)
  
  ## Duration calc
  f <- "odin/ode_anp_cohort.R"
  model <- odin::odin(f, target = "r")
  cm <- model()
  calc_dur <- function(cm, inp) {
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
      All = a + sn + sp
    )
  }
  
  
  res <- lapply(c("female", "male"), function(sex) {
    dat <- c(dat_anp[[sex]], inp)
    fitted <- sampling(m_at_ts, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain, data = dat)
    mc <- extract(fitted, pars = c("del_sn", "del_sp", "r_tr", "adr", "p_sp"))

    duration <- apply(as.data.frame(mc)[c("del_sn", "del_sp", "r_tr", "p_sp")], 1, function(x) {
      calc_dur(cm, c(list(r_sym = 12 / dur_asym, A0 = 1, Sn0 = 0, Sp0 = 0), x))
    })
    duration <- t(duration)
    
    list(
      DurAsym = dur_asym,
      MC = mc,
      Parameters = summary(fitted, pars = c("del_sn", "del_sp", "r_tr", "adr", "p_sp"))$summary,
      Proportion = summary(fitted, pars = c("pr_a", "pr_sn", "pr_sp"))$summary,
      Prevalence = summary(fitted, pars = "prv")$summary,
      Incidence = list(
        Asym = summary(fitted, pars = "inc_a")$summary,
        Sym = summary(fitted, pars = "inc_s")$summary,
        Sn = summary(fitted, pars = "inc_sn")$summary,
        Sp = summary(fitted, pars = "inc_sp")$summary
      ),
      Duration = duration,
      Data = c(
        PrvSn = (dat$Sn / dat$N),
        PrvSp = (dat$Sp / dat$N),
        NotiSn = (dat$Noti_Sn / dat$Pop)[1],
        NotiSp = (dat$Noti_Sp / dat$Pop)[1],
        DelaySn = (dat$Sn / dat$N) / (dat$Noti_Sn / dat$Pop)[1],
        DelaySp = (dat$Sp / dat$N) / (dat$Noti_Sp / dat$Pop)[1]
      ),
      Exo = inp
    )
  })
  names(res) <- c("female", "male")
  res
}


fit_np <- function(dat_np, n_iter = 5E4, n_collect = 1E3, n_chain = 3) {
  require(rstan)
  
  m_t_ts <- readRDS(file = "stan/m_t_ts.rds")
  
  pars <- c("del_sp", "del_sn", "r_tr", "adr", "p_sp", "pr_sn", "pr_sp",
            "inc_s", "inc_sn", "inc_sp", "noti_sn", "noti_sp", "prv")
  
  
  inp <- list(r_death_sp = 0.127, r_death_sn = 0.024, r_sc = 0.15)
  
  ## Duration calc
  f <- "odin/ode_np_cohort.R"
  model <- odin::odin(f, target = "r")
  cm <- model()
  calc_dur <- function(cm, inp) {
    cm$set_user(user = inp)
    res <- cm$run(seq(0, 100, 0.1))
    
    sn <- sum((res[-1, 2] + res[-nrow(res), 2]) * diff(res[, "t"]) / 2)
    sp <- sum((res[-1, 3] + res[-nrow(res), 3]) * diff(res[, "t"]) / 2)
    c(
      A = 0,
      Sn = sn,
      Sp = sp,
      Sym = sn + sp,
      All = sn + sp
    )
  }
  
  
  res <- lapply(c("female", "male"), function(sex) {
    dat <- c(dat_np[[sex]], inp)
    fitted <- sampling(m_t_ts, iter = n_iter, warmup = n_iter - n_collect, chain = n_chain, data = dat)
    
    mc <- extract(fitted, pars = c("del_sn", "del_sp", "r_tr", "adr", "p_sp"))
    
    duration <- apply(as.data.frame(mc)[c("del_sn", "del_sp", "r_tr", "p_sp")], 1, function(x) {
      calc_dur(cm, c(list(n0  = 1), x))
    })
    duration <- t(duration)
    
    list(
      DurAsym = 0,
      MC = mc,
      Parameters = summary(fitted, pars = c("del_sn", "del_sp", "r_tr", "adr", "p_sp"))$summary,
      Proportion = summary(fitted, pars = c("pr_sn", "pr_sp"))$summary,
      Prevalence = summary(fitted, pars = "prv")$summary,
      Incidence = list(
        Sym = summary(fitted, pars = "inc_s")$summary,
        Sn = summary(fitted, pars = "inc_sn")$summary,
        Sp = summary(fitted, pars = "inc_sp")$summary
      ),
      Duration = duration,
      Data = c(
        PrvSn = (dat$Sn / dat$N),
        PrvSp = (dat$Sp / dat$N),
        NotiSn = (dat$Noti_Sn / dat$Pop)[1],
        NotiSp = (dat$Noti_Sp / dat$Pop)[1],
        DelaySn = (dat$Sn / dat$N) / (dat$Noti_Sn / dat$Pop)[1],
        DelaySp = (dat$Sp / dat$N) / (dat$Noti_Sp / dat$Pop)[1]
      ),
      Exo = inp
    )
  })
  names(res) <- c("female", "male")
  res
}
