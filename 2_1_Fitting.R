rm(list = ls())


source("R/fit.R")

library(rstan)
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = "-march=corei7 -mtune=corei7")



n_iter <- 1E4
n_collect <- 500
n_chain <- 2
prior_anp <- list(r_death_sp = 0.12, r_death_sn = 0.022, r_sc_l = 0.1, r_sc_u = 0.3)
prior_as <- list(r_death = 0.071, r_sc_l = 0.1, r_sc_u = 0.3)
prior_np <- list(r_death_sp = 0.12, r_death_sn = 0.022, r_sc_l = 0.1, r_sc_u = 0.3)
prior_asc <- prior_as <- list(r_death = 0.071, r_sc_l = 0.1, r_sc_u = 0.3)


load("data/Input_Kenya.rdata")


d <- fit_anp(prv$a_sn_sp$Female, n_iter = n_iter, n_chain = n_chain, prior = prior_anp)

d <- fit_as(prv$a_sym$Male, n_iter = n_iter, n_chain = n_chain, prior = prior_as)

d <- fit_np(prv$sn_sp$Male, n_iter = n_iter, n_chain = n_chain, prior = prior_np)

d <- fit_asc(prv$a_sym_cs$T_25_34, n_iter = n_iter, n_chain = n_chain, prior = prior_asc)
