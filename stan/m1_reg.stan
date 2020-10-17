data {
  // Data from a prevalence survey
  int n_gp;
  int<lower=0> N[n_gp];
  int<lower=0> Asym[n_gp];
  int<lower=0> Sn[n_gp];
  int<lower=0> Sp[n_gp];
  real YearSurveyed; // timing of the survey


  // Time-series of notification data
  int<lower=0> n_t; // number of time points
  int<lower=0> Pop[n_gp, n_t]; // population size
  int<lower=0> NotiSn[n_gp, n_t]; // sm- counts
  int<lower=0> NotiSp[n_gp, n_t]; // sm+ counts
  real<lower=0> Years[n_t]; // years of the notification data

  // Prior knowledge
  real<lower=0> r_sc_l;
  real<lower=0> r_sc_u;
  real<lower=0> scale_dur;

  // Exogenous variables
  real<lower=0> r_death_bg[n_gp];
  real<lower=0> r_death_sn[n_gp];
  real<lower=0> r_death_sp[n_gp];
}
parameters {
  real<lower=0, upper=1> p_sp; // Probability of sm+ at symptom onset
  real<lower=0> r_tr; // conversion rate, sm- to sm+
  real<lower=r_sc_l, upper = r_sc_u> r_sc; // self-cure rate
  
  real<lower=0> r_det_sn0;
  real<lower=0> r_det_sp0;
  real<lower=0> r_sym0;
  
  vector[n_gp] lrr_cs_sp;
  vector[n_gp] lrr_cs_sn;
  vector[n_gp] lrr_sym;
  
  vector<lower=0, upper=1>[n_gp] prv0; // prevalence of all active TB
  vector<lower=-0.15, upper=0.15>[n_gp] adr; // annual decline rate
}
transformed parameters {
  vector<lower=0>[n_gp] r_det_sn; // delay to notification, sm-
  vector<lower=0>[n_gp] r_det_sp; // delay to notification, sm+
  vector<lower=0>[n_gp] r_sym; // symptom development rate == inversion of asymptomatic phase
  
  
  vector<lower=0>[n_gp] ra;
  vector<lower=0>[n_gp] rn;
  vector<lower=0>[n_gp] rp;

  vector<lower=0>[n_gp] sn0;
  vector<lower=0>[n_gp] sp0;
  vector<lower=0, upper=1>[n_gp] pr_a;
  vector<lower=0, upper=1>[n_gp] pr_sn;
  vector<lower=0, upper=1>[n_gp] pr_sp;
  vector<lower=0, upper=1>[n_gp] prv_a;
  vector<lower=0, upper=1>[n_gp] prv_sn;
  vector<lower=0, upper=1>[n_gp] prv_sp;

  matrix<lower=0>[n_gp, n_t] prv;
  matrix<lower=0>[n_gp, n_t] nr_sn;
  matrix<lower=0>[n_gp, n_t] nr_sp;
  matrix<lower=0>[n_gp, n_t] nr;
  
  r_det_sn[1] = r_det_sn0;
  r_det_sp[1] = r_det_sp0;
  r_sym[1] = r_sym0;
  
  for (i in 2:n_gp) {
    r_det_sn[i] = r_det_sn0 * exp(lrr_cs_sn[i]);
    r_det_sp[i] = r_det_sp0 * exp(lrr_cs_sp[i]);
    r_sym[i] = r_sym0 * exp(lrr_sym[i]);
  }
  
  

  for (j in 1:n_gp) {
    ra[j] = r_sc + r_death_bg[j];
    rn[j] = r_sc + r_death_bg[j] + r_death_sn[j];
    rp[j] = r_sc + r_death_bg[j] + r_death_sp[j];
    
    sn0[j] = (1 - p_sp) * r_sym[j] / (r_tr + rn[j] + r_det_sn[j] - adr[j]);
    sp0[j] = (p_sp * r_sym[j] + r_tr * sn0[j]) / (rp[j] + r_det_sp[j] - adr[j]);

    // probabilities of the states among all active TB
    pr_a[j] = 1 / (1 + sn0[j] + sp0[j]);
    pr_sn[j] = sn0[j] * pr_a[j];
    pr_sp[j] = sp0[j] * pr_a[j];

    // prevalence of the states
    prv_a[j] = prv0[j] * pr_a[j];
    prv_sn[j] = prv0[j] * pr_sn[j];
    prv_sp[j] = prv0[j] * pr_sp[j];
  }


  for (i in 1:n_t) {
    for (j in 1:n_gp) {
      // forecasts/backcasts of prevalence
      prv[j, i] = prv0[j] * exp(- adr[j] * (Years[i] - YearSurveyed));

      // notification rates
      nr_sn[j, i] = prv[j, i] * pr_sn[j] * r_det_sn[j];
      nr_sp[j, i] = prv[j, i] * pr_sp[j] * r_det_sp[j];
      nr[j, i] = nr_sn[j, i] + nr_sp[j, i];
    }
  }
}
model {
  p_sp ~ uniform(0, 1);
  r_tr ~ uniform(0, 0.5);
  r_sc ~ uniform(r_sc_l, r_sc_u);

  prv0 ~ uniform(0, 1);
  r_det_sp0 ~ inv_gamma(scale_dur, scale_dur);
  r_det_sn0 ~ inv_gamma(scale_dur, scale_dur);
  r_sym0 ~ inv_gamma(scale_dur, scale_dur);
  adr ~ uniform(-0.15, 0.15);
  
  for (i in 1:n_gp) {
    lrr_cs_sp[i] ~ normal(0, 1);
    lrr_cs_sn[i] ~ normal(0, 1);
    lrr_sym[i] ~ normal(0, 1);
  }
  
  for (j in 1:n_gp) {
    // prevalence to prevalence survey data
    target += binomial_lpmf(Asym[j] | N[j], prv_a[j]);
    target += binomial_lpmf(Sn[j] | N[j], prv_sn[j]);
    target += binomial_lpmf(Sp[j] | N[j], prv_sp[j]);

    // notification rate to notification data
    for (i in 1:n_t) {
      target += poisson_lpmf(NotiSn[j, i] | nr_sn[j, i] * Pop[j, i]);
      target += poisson_lpmf(NotiSp[j, i] | nr_sp[j, i] * Pop[j, i]);
    }
  }
}
generated quantities {
  real log_lik_noti[n_t, 2, n_gp];
  real log_lik_pr[3, n_gp];

  for (j in 1:n_gp) {
    // prevalence to prevalence survey data
    log_lik_pr[1, j] = binomial_lpmf(Asym[j] | N[j], prv_a[j]);
    log_lik_pr[2, j] = binomial_lpmf(Sn[j] | N[j], prv_sn[j]);
    log_lik_pr[3, j] = binomial_lpmf(Sp[j] | N[j], prv_sp[j]);

    // notification rate to notification data
    for (i in 1:n_t) {
      log_lik_noti[i, 1, j] = poisson_lpmf(NotiSn[j, i] | nr_sn[j, i] * Pop[j, i]);
      log_lik_noti[i, 2, j] = poisson_lpmf(NotiSp[j, i] | nr_sp[j, i] * Pop[j, i]);
    }
  }
}
