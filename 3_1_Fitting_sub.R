rm(list = ls())


source("R/fit.R")

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores() - 2)

n_iter <- 1E4
n_collect <- 500
n_chain <- 2
prior_anp <- list(r_death_sp = 0.12, r_death_sn = 0.022, r_sc_l = 0.1, r_sc_u = 0.3)
prior_as <- list(r_death = 0.071, r_sc_l = 0.1, r_sc_u = 0.3)
prior_np <- list(r_death_sp = 0.12, r_death_sn = 0.022, r_sc_l = 0.1, r_sc_u = 0.3)
prior_asc <- prior_as <- list(r_death = 0.071, r_sc_l = 0.1, r_sc_u = 0.3)


countries <- c(
  "Cambodia", "Lao People's Democratic Republic", 
  "Malawi", "Pakistan", "Philippines", 
  "United Republic of Tanzania", "Uganda" , "Viet Nam", "Zambia"
)

countries_cs <- c(
  "Malawi", "Philippines", 
  "United Republic of Tanzania", "Uganda" , "Zambia"
)



for (country in countries) {
  load(paste0("data/Input_", country, ".rdata"))
  
  post_a_sn_sp <- lapply(prv$a_sn_sp, fit_anp, n_iter = n_iter, n_chain = n_chain, prior = prior_anp)
  post_a_sym <- lapply(prv$a_sym, fit_as, n_iter = n_iter, n_chain = n_chain, prior = prior_as)
  #post_sn_sp <- lapply(prv$sn_sp, fit_np, n_iter = n_iter, n_chain = n_chain, prior = prior_np)
  
  if (country %in% countries_cs) {
    post_a_sym_cs <- lapply(prv$a_sym_cs, fit_asc, n_iter = n_iter, n_chain = n_chain, prior = prior_asc)
    save(post_a_sym, post_a_sym_cs, file =paste0("output/Post_", country, ".rdata"))
  } else {
    save(post_a_sym, file =paste0("output/Post_", country, ".rdata"))
  }
  
}


