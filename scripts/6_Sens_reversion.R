library(tidyverse)


test_est <- function(n_sim = 500, r_rev = 0, r_seek = 0.1, r_sym = 0.1, pars = list()) {
  model <- odin::odin({
    deriv(A) <- inc - r_sym * A - ga_a * A + r_rev * S
    deriv(S) <- r_sym * A - r_rev * S - r_seek * S - ga_s * S
    
    
    initial(A) <- A0 * N * prv0
    initial(S) <- S0 * N * prv0
    
    output(PrvA) <- A / N
    output(PrvS) <- S / N
    output(CNR) <- r_seek * S / N  
    
    ga_a <- r_sc + r_death
    ga_s <- r_sc + r_death + r_death_s
    
    inc <- (ga_a * A + ga_s * S + r_seek * S) - adr * (A + S)
    
    
    A0 <- (r_seek + ga_s + r_rev - adr) / (r_sym + r_seek + ga_s + r_rev - adr)
    S0 <- r_sym / (r_sym + r_seek + ga_s + r_rev - adr)
    
    N <- 1E7
    prv0 <- 0.01
    
    # Parameters
    adr <- user(0.01)
    r_sc <- user(0.2) # Self-cure
    r_death <- user(0.015)
    r_death_s <- user(0.1)
    
    r_rev <- user(0)
    r_sym <- user(1) # Symptom activation
    r_seek <- user(1) # Delay sp 
  })
  
  
  fn_post <- function(x, dat, pars) {
    r_sym <- x[1]
    r_seek <- x[2]
    prv <- x[3]
    
    li <- with(c(dat, pars), {
      ga_a <- r_sc + r_death
      ga_s <- ga_a + r_death_s
      
      asym <- (r_seek + ga_s - adr) / (r_sym + r_seek + ga_s - adr) * prv
      sym <- r_sym / (r_sym + r_seek + ga_s - adr) * prv
      
      li <- dgamma(1 / r_sym, 1, 1, log = T)
      li <- li + dgamma(1 / r_seek, 1, 1, log = T)
      li <- li + dbinom(Asym, Pop, asym, log = T)
      li <- li + dbinom(Sym, Pop, sym, log = T)
      
      li <- li + sum(dpois(Noti, Pop * r_seek * sym * exp(- adr * (Year - SurveyYear)), log = T))
      li
    })
    
    - li
  }
  
  
  cm <- model(user = pars)
  cm$set_user(r_rev = r_rev)
  cm$set_user(r_sym = r_sym)
  cm$set_user(r_seek = r_seek)
  
  ys <- cm$run(seq(0, 10, .1))
  ys <- as.data.frame(ys[ys[, 1] == round(ys[, 1]), ])[-1, ]
  
  est <- bind_rows(lapply(1:500, function(i) {
    dat <- list(
      Year = 1:10,
      SurveyYear = 5,
      Pop = 1E7,
      Asym = rbinom(1, ys$PrvA[ys$t == 5], size = 1E7),
      Sym = rbinom(1, ys$PrvS[ys$t == 5], size = 1E7),
      Noti = rpois(10, ys$CNR * 1E7)
    )
    
    fit <- nlminb(c(0.01, 0.01, 0.001), fn_post, 
                  lower = 1E-6, upper = c(10, 10, 1), 
                  dat = dat, pars = pars)
    x <- fit$par
    x <- c(i, x)
    names(x) <- c("Key", "r_sym", "r_seek", "prv")
    as.list(x)
  }))
  
  
  est <- est %>% 
    mutate(
      DurA = 1 / (r_sym + pars$r_sc + pars$r_death) ,
      DurS = 1 / (r_seek + pars$r_sc + pars$r_death + pars$r_death_s)
    )
  
  est$r_sym_true <- r_sym
  est$r_seek_true <- r_seek
  est$r_rev_true <- r_rev
  
  est$DurA_true <- 1 / (r_sym + pars$r_sc + pars$r_death)
  est$DurS_true <- 1 / (r_seek + r_rev + pars$r_sc + pars$r_death + pars$r_death_s)
  est
} 



pars <- list(
  adr = 0.01,
  r_sc = 0.2,
  r_death = 0.015,
  r_death_s = 0.1
)


est <- bind_rows(lapply(c(0, 0.1, 0.2, 0.3), function(r_rev) {
  bind_rows(lapply(seq(0.5, 1.5, 0.5), function(r_sym) {
    bind_rows(lapply(seq(0.5, 1.5, 0.5), function(r_seek) {
      test_est(500, r_rev = r_rev, r_sym = r_sym, r_seek = r_seek, pars = pars)
    }))
  }))
}))


stat <- est %>% 
  mutate(err_r_sym = r_sym / r_sym_true - 1, err_r_seek = r_seek / r_seek_true - 1,
         err_DurA = DurA / DurA_true - 1, err_DurS = DurS / DurS_true - 1) %>% 
  group_by(r_rev_true, r_sym_true, r_seek_true, DurA_true, DurS_true) %>% 
  summarise(across(-Key, list(M = mean, 
                              L = function(x) quantile(x, 0.025),
                              U = function(x) quantile(x, 0.975)))) %>% 
  pivot_longer(-c(r_rev_true, r_sym_true, r_seek_true)) %>% 
  tidyr::extract(name, c("Index", "Stat"), "(\\w+)_(M|L|U)") %>% 
  mutate(
    r_sym_exp = as.factor(paste("theta", r_sym_true, sep = ": ")),
    r_seek_exp = as.factor(paste("rho", r_seek_true, sep = ": ")),
  )
  


save(est, stat, file = "out/Sens/Reversion.rdata")



## Visualisation

library(ggplot2)

theme_set(theme_bw() + theme(text = element_text(family = "serif")))

gs <- list()


gs$g_est <- stat %>% 
  filter(Index %in% c("r_sym", "r_seek")) %>% 
  pivot_wider(values_from = value, names_from = Stat) %>% 
  ggplot(aes(x = r_rev_true)) +
  geom_line(aes(y = M, colour = Index)) +
  geom_ribbon(aes(ymin = L, ymax = U, fill = Index), alpha = 0.3) + 
  scale_y_continuous("MAP estimator (per year)", breaks = seq(0, 1.5, 0.5)) +
  scale_x_continuous("Symptom reversion rate, per year") +
  scale_colour_discrete("Rate estimator", labels = c(r_sym = "Symptom development", r_seek = "Care-seeking")) +
  facet_wrap(r_sym_exp + r_seek_exp~., nrow = 3, labeller = label_parsed) +
  guides(fill = "none") +
  expand_limits(y = 0) +
  theme(legend.position = "bottom")


gs$g_err_rate <- stat %>% 
  filter(Index %in% c("err_r_sym", "err_r_seek")) %>% 
  pivot_wider(values_from = value, names_from = Stat)%>% 
  ggplot(aes(x = r_rev_true)) +
  geom_line(aes(y = M, colour = Index)) +
  geom_ribbon(aes(ymin = L, ymax = U, fill = Index), alpha = 0.3) + 
  scale_y_continuous("Percentage error, (Estimeted - Actual) / Actual * 100%", label = scales::percent_format(1)) +
  scale_x_continuous("Symptom reversion rate, per year") +
  scale_color_discrete("Rate estimator", labels = c(err_r_sym = "Symptom development", err_r_seek = "Care-seeking")) +
  facet_wrap(r_sym_exp + r_seek_exp~., nrow = 3, labeller = label_parsed) +
  guides(fill = "none") +
  expand_limits(y = - .5) +
  theme(legend.position = "bottom")


gs$g_err_dur <- stat %>% 
  filter(Index %in% c("err_DurS", "err_DurA")) %>% 
  pivot_wider(values_from = value, names_from = Stat)%>% 
  ggplot(aes(x = r_rev_true)) +
  geom_line(aes(y = M, colour = Index)) +
  geom_ribbon(aes(ymin = L, ymax = U, fill = Index), alpha = 0.3) + 
  scale_y_continuous("Prediction error, (Estimeted - Actual) / Actual * 100%", label = scales::percent_format(1)) +
  scale_x_continuous("Symptom reversion rate, per year") +
  scale_color_discrete("Duration estimator", labels = c(err_DurS = "Symptomatic stage", err_DurA = "Asymptomatic stage")) +
  facet_wrap(r_sym_exp + r_seek_exp~., nrow = 3, labeller = label_parsed) +
  guides(fill = "none") +
  expand_limits(y = .5) +
  theme(legend.position = "bottom")


out_path <- glue::as_glue("docs/figs/sens/")

ggsave(gs$g_est, filename = out_path + "g_est.pdf", width = 7, height = 8)
ggsave(gs$g_err_rate, filename = out_path + "g_err_rate.pdf", width = 7, height = 8)
ggsave(gs$g_err_dur, filename = out_path + "g_err_dur.pdf", width = 7, height = 8)



