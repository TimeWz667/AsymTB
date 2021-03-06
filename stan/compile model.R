## Compile models -----
library(rstan)
rstan_options(auto_write = TRUE)


m0 <- stan_model("stan/m0.stan")


# Marginalised
m1_uni <- stan_model("stan/m1_uni.stan")

m1_dur <- stan_model("stan/m1_duration_free.stan")

m1_fixed <- stan_model("stan/m1_fixed.stan")

m1_reg <- stan_model("stan/m1_reg.stan")


# Multi regression
m2_cov <- stan_model("stan/m2_cov.stan")


# 
m3_as_uni <- stan_model("stan/m3_as_uni.stan")
m3_as_fixed <- stan_model("stan/m3_as_fixed.stan")
m3_as_dur <- stan_model("stan/m3_as_duration_free.stan")

m3_asc_uni <- stan_model("stan/m3_asc_uni.stan")
m3_asc_fixed <- stan_model("stan/m3_asc_fixed.stan")
m3_asc_dur <- stan_model("stan/m3_asc_duration_free.stan")
# m3_anp <- stan_model("stan/m3_anp.stan")
