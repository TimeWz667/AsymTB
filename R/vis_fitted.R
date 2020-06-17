
visualise_fitted <- function(fit_asym, fit_baseline, inc) {
  require("ggplot2")
  
  ds <- list()
  
  fitted <- c(
    list(fit_baseline),
    fit_asym
  )
  
  rat_pn <- rbind(
    data.frame(Smear = c("Smear-", "Smear+"), Sex = "Female", PrvNoti = fit_baseline$female$Data[5:6]),
    data.frame(Smear = c("Smear-", "Smear+"), Sex = "Male", PrvNoti = fit_baseline$male$Data[5:6])
  ) 
  
  rat_apn <- rbind(
    data.frame(Smear = c("Smear-", "Smear+"), Sex = "Female", PrvNoti = fit_asym[[1]]$female$Data[5:6]),
    data.frame(Smear = c("Smear-", "Smear+"), Sex = "Male", PrvNoti = fit_asym[[1]]$male$Data[5:6])
  )
  
  dat <- data.table::rbindlist(lapply(fitted, function(x) {
    d <- rbind(data.frame(s=c("Smear-", "Smear+"), Sex = "Female", x$female$Parameters[c("del_sn", "del_sp"), c(4:8)]),
               data.frame(s=c("Smear-", "Smear+"), Sex = "Male", x$male$Parameters[c("del_sn", "del_sp"), c(4:8)])
    )
    colnames(d) <- c("Smear", "Sex", "l95", "l50", "m", "u50", "u95")
    d$Duration <- x$female$DurAsym
    d
  }))
  
  
  g_delay <- ggplot(dat, aes(x = Duration)) +
    geom_crossbar(aes(y = m, ymin = l95, ymax = u95, fill = "95%"), width = 1, colour = "white") +
    geom_crossbar(aes(y = m, ymin = l50, ymax = u50, fill = "50%"), width = 1, colour = "white") +
    geom_hline(data = rat_pn, aes(yintercept = PrvNoti, linetype = "Noti / Prev")) +
    geom_hline(data = rat_apn, aes(yintercept = PrvNoti, linetype = "Noti / Prev, Sym")) +
    scale_fill_discrete("Quantiles") +
    scale_x_continuous("Asymptomatic phase (months)", breaks = c(0, 3, 6, 9, 12)) +
    scale_y_continuous("Delay to treatment from symptom onset (years)") +
    facet_grid(Smear ~ Sex) +
    expand_limits(y = c(0, NA)) +
    theme(legend.position = "bottom")
  
  
  dat <- data.table::rbindlist(lapply(fitted, function(x) {
    d <- rbind(data.frame(s=c("Smear-", "Smear+"), Sex = "Female", 
                          x$female$Parameters[c("del_sn", "del_sp"), c(4:8)] + x$female$DurAsym / 12),
               data.frame(s=c("Smear-", "Smear+"), Sex = "Male", 
                          x$male$Parameters[c("del_sn", "del_sp"), c(4:8)] + x$female$DurAsym / 12)
    )
    colnames(d) <- c("Smear", "Sex", "l95", "l50", "m", "u50", "u95")
    d$Duration <- x$female$DurAsym
    d
  }))
  
  g_err <- ggplot(dat, aes(x = Duration)) +
    geom_crossbar(aes(y = m, ymin = l95, ymax = u95, fill = "95%"), width = 1, colour = "white") +
    geom_crossbar(aes(y = m, ymin = l50, ymax = u50, fill = "50%"), width = 1, colour = "white") +
    geom_hline(data = rat_pn, aes(yintercept = PrvNoti, linetype = "Noti / Prev")) +
    geom_hline(data = rat_apn, aes(yintercept = PrvNoti, linetype = "Noti / Prev, Sym")) +
    scale_fill_discrete("Quantiles") +
    scale_linetype_discrete("") +
    scale_x_continuous("Asymptomatic phase (months)", breaks = c(0, 3, 6, 9, 12)) +
    scale_y_continuous("Asymptomatic phase + Care-seeking delay (years)") +
    facet_grid(Smear ~ Sex) +
    expand_limits(y = c(0, NA)) +
    theme(legend.position = "bottom")
  
  
  dat <- data.table::rbindlist(lapply(fitted, function(x) {
    d <- rbind(
      quantile(x$female$Duration[, "All"], c(0.025, 0.25, 0.5, 0.75, 0.975)),
      quantile(x$male$Duration[, "All"], c(0.025, 0.25, 0.5, 0.75, 0.975))
    )
    colnames(d) <- c("l95", "l50", "m", "u50", "u95")
    
    d <- data.frame(Sex = c("Female", "Male"), d)
    d$Duration <- x$female$DurAsym
    d
  }))
  
  ds$Duration <- dat
  
  g_dur <- ggplot(dat, aes(x = Duration)) +
    geom_crossbar(aes(y = m, ymin = l95, ymax = u95, fill = "95%"), width = 1, colour = "white") +
    geom_crossbar(aes(y = m, ymin = l50, ymax = u50, fill = "50%"), width = 1, colour = "white") +
    #geom_hline(data = rat_pn, aes(yintercept = PrvNoti, linetype = "Noti / Prev")) +
    #geom_hline(data = rat_apn, aes(yintercept = PrvNoti, linetype = "Noti / Prev, Sym")) +
    scale_fill_discrete("Quantiles") +
    scale_linetype_discrete("") +
    scale_x_continuous("Asymptomatic phase, months", breaks = c(0, 3, 6, 9, 12)) +
    scale_y_continuous("Duration of pre-diagnostic phase, years") +
    facet_grid(. ~ Sex) +
    expand_limits(y = c(0, NA)) +
    theme(legend.position = "bottom")
  
  
  dat <- data.table::rbindlist(lapply(fitted, function(x) {
    d <- rbind(data.frame(Sex = "Female", Inc = "Symptomatic", t(x$female$Incidence$Sym[3, c(4:8)])),
               data.frame(Sex = "Male", Inc = "Symptomatic", t(x$male$Incidence$Sym[3, c(4:8)]))
    )
    
    if (length(x$female$Incidence$Asym)) {
      d <- rbind(
        d,
        data.frame(Sex = "Female", Inc = "Asymptomatic", t(x$female$Incidence$Asym[3, c(4:8)])),
        data.frame(Sex = "Male", Inc = "Asymptomatic", t(x$male$Incidence$Asym[3, c(4:8)]))
      )
    }
    
    colnames(d) <- c("Sex", "Inc", "l95", "l50", "m", "u50", "u95")
    d$Duration <- x$female$DurAsym
    d
  }))
  
  ds$Inc <- dat
  ds$IncD <- inc
  
  
  g_inc <- ggplot(dat, aes(x = Duration)) +
    geom_pointrange(aes(y = m, ymin = l95, ymax = u95, colour = Inc), 
                    position = position_dodge(width = 0.5), size = rel(0.5)) +
    geom_pointrange(data = inc, aes(x = -1, y = m, ymin = l, ymax = u, colour = "WHO estimates"), 
                    position = position_dodge(width = 0.5)) +
    scale_color_discrete("Quantiles") +
    scale_linetype_discrete("") +
    scale_x_continuous("Asymptomatic phase, months", breaks = c(0, 3, 6, 9, 12)) +
    scale_y_continuous("Incidence rate, per 100 000", labels = function(x) x * 1E5) +
    facet_wrap(Sex ~ .) +
    expand_limits(y = 0) +
    theme(legend.position = "bottom")
  
  list(
    delay = g_delay,
    error = g_err,
    duration = g_dur,
    inc = g_inc,
    dat = ds
  )
}


visualise_bind <- function(gf) {
  duration <- data.table::rbindlist(lapply(names(gf), function(country) {
    dat <- gf[[country]]$dat$Duration
    dat$GP = country
    dat
  }))
  
  g_dur <- ggplot(duration, aes(x = Duration)) +
    geom_crossbar(aes(y = m, ymin = l95, ymax = u95, fill = "95%"), width = 1, colour = "white") +
    geom_crossbar(aes(y = m, ymin = l50, ymax = u50, fill = "50%"), width = 1, colour = "white") +
    #geom_hline(data = rat_pn, aes(yintercept = PrvNoti, linetype = "Noti / Prev")) +
    #geom_hline(data = rat_apn, aes(yintercept = PrvNoti, linetype = "Noti / Prev, Sym")) +
    geom_vline(xintercept = 0.5) +
    scale_fill_discrete("Quantiles") +
    scale_linetype_discrete("") +
    scale_x_continuous("Asymptomatic phase, months", breaks = c(0, 3, 6, 9, 12)) +
    scale_y_continuous("Duration of pre-diagnostic phase, years") +
    facet_grid(GP ~ Sex) +
    expand_limits(y = c(0, NA)) +
    theme(legend.position = "bottom")
  
  inc <- data.table::rbindlist(lapply(names(gf), function(country) {
    dat <- gf[[country]]$dat$Inc
    dat$GP = country
    dat
  }))
  
  incd <- data.table::rbindlist(lapply(names(gf), function(country) {
    dat <- gf[[country]]$dat$IncD
    dat$GP = country
    dat
  }))
  
  
  g_inc <- ggplot(inc, aes(x = Duration)) +
    geom_pointrange(aes(y = m, ymin = l95, ymax = u95, colour = Inc), 
                    position = position_dodge(width = 0.5), size = rel(0.5)) +
    geom_pointrange(data = incd, aes(x = -1, y = m, ymin = l, ymax = u, colour = "WHO estimates"), 
                    position = position_dodge(width = 0.5)) +
    scale_color_discrete("Incidences") +
    scale_linetype_discrete("") +
    scale_x_continuous("Asymptomatic phase, months", breaks = c(0, 3, 6, 9, 12)) +
    scale_y_continuous("Incidence rate, per 100 000", labels = function(x) x * 1E5) +
    facet_grid(GP ~ Sex) +
    expand_limits(y = 0) +
    theme(legend.position = "bottom")
  
  list(duration = g_dur,
       inc = g_inc)
}

