data {
  // Data from a prevalence survey
  int n_gp;
  int<lower=0> N[n_gp];
  int<lower=0> Asym[n_gp];
  int<lower=0> Sym[n_gp];
  real YearSurveyed; // timing of the survey


  // Time-series of notification data
  int<lower=0> n_t; // number of time points
  int<lower=0> Pop[n_gp, n_t]; // population size
  int<lower=0> Noti[n_gp, n_t]; // notification counts
  real<lower=0> Years[n_t]; // years of the notification data

  // Prior knowledge
  real<lower=0> r_sc_l;
  real<lower=0> r_sc_u;
  real<lower=0> scale_dur;
  real<lower=0, upper=1> pr_fp;

  // Exogenous variables
  real<lower=0> r_death_a[n_gp];
  real<lower=0> r_death_s[n_gp];
}
parameters {
  vector<lower=0>[n_gp] r_det;
  vector<lower=0, upper=1>[n_gp] prv0;
  vector<lower=-0.2, upper=0.2>[n_gp] adr;
  vector<lower=0>[n_gp] r_sym;
  real<lower=r_sc_l, upper = r_sc_u> r_sc;
}
transformed parameters {
  vector<lower=0>[n_gp] ra;
  vector<lower=0>[n_gp] rs;

  vector<lower=0>[n_gp] a0;
  vector<lower=0>[n_gp] s0;
  
  vector<lower=0, upper=1>[n_gp] pr_a;
  vector<lower=0, upper=1>[n_gp] pr_s;
  vector<lower=0, upper=1>[n_gp] prv_a;
  vector<lower=0, upper=1>[n_gp] prv_s;

  matrix<lower=0>[n_gp, n_t] prv;
  matrix<lower=0>[n_gp, n_t] nr;
  
  
  for (j in 1:n_gp) {
    ra[j] = r_sc + r_death_a[j];
    rs[j] = r_sc + r_death_s[j];
    
    a0[j] = rs[j] + r_det[j] - adr[j];
    s0[j] = r_sym[j];

    // probabilities of the states among all active TB
    pr_a[j] = a0[j] / (a0[j] + s0[j]);
    pr_s[j] = s0[j] / (a0[j] + s0[j]);

    // prevalence of the states
    prv_a[j] = prv0[j] * pr_a[j];
    prv_s[j] = prv0[j] * pr_s[j];
  }
  
  for (i in 1:n_t) {
    for (j in 1:n_gp) {
      prv[j, i] = prv0[j] * exp(- adr[j] * (Years[i] - YearSurveyed));
      nr[j, i] = prv[j, i] * pr_s[j] * r_det[j];
    }
  }
}
model {
  prv0 ~ uniform(0, 1);
  r_det ~ inv_gamma(scale_dur, scale_dur);

  r_sym ~ inv_gamma(scale_dur, scale_dur);
  r_sc ~ uniform(r_sc_l, r_sc_u);
  
  adr ~ uniform(-0.2, 0.2);

  for (j in 1:n_gp) {
    target += binomial_lpmf(Asym[j] | N[j], prv_a[j]);
    target += binomial_lpmf(Sym[j] | N[j], prv_s[j]);
    
    for (i in 1:n_t) {
      target += poisson_lpmf(Noti[j, i] | nr[j, i] * Pop[j, i] / (1 - pr_fp));
    }
  }
}
generated quantities {
  real log_lik_noti[n_t, n_gp];
  real log_lik_pr[2, n_gp];
  
  real<lower=0> IncN_A[n_gp];
  real<lower=0> IncN_S[n_gp];
  real<lower=0> NotiN[n_gp];
  real<lower=0> PrvN_A[n_gp];
  real<lower=0> PrvN_S[n_gp];
  real<lower=0> CDR_A[n_gp];
  real<lower=0> CDR_S[n_gp];
  real<lower=0> Gap_A[n_gp];
  real<lower=0> Gap_S[n_gp];
  
  real dur_a[n_gp];
  real dur_s[n_gp];

  for (j in 1:n_gp) {
    // prevalence to prevalence survey data
    log_lik_pr[1, j] = binomial_lpmf(Asym[j] | N[j], prv_a[j]);
    log_lik_pr[2, j] = binomial_lpmf(Sym[j] | N[j], prv_s[j]);

    // notification rate to notification data
    for (i in 1:n_t) {
      log_lik_noti[i, j] = poisson_lpmf(Noti[j, i] | nr[j, i] * Pop[j, i]);
    }
    
    
    NotiN[j] = nr[j, n_t] * Pop[j, n_t];
    
    PrvN_A[j] = pr_a[j] * prv[j, n_t] * Pop[j, n_t];
    PrvN_S[j] = pr_s[j] * prv[j, n_t] * Pop[j, n_t];
    
    IncN_A[j] = (ra[j] * pr_a[j] + (rs[j] + r_det[j]) * pr_s[j] - adr[j]) * prv[j, n_t] * Pop[j, n_t];
    IncN_S[j] = r_sym[j] * PrvN_A[j];
    
    CDR_A[j] = NotiN[j] / IncN_A[j];
    CDR_S[j] = NotiN[j] / IncN_S[j];
    
    dur_a[j] = 1 / (ra[j] + r_sym[j]);
    dur_s[j] = 1 / (rs[j] + r_det[j]);
    Gap_A[j] = r_sym[j] * dur_a[j];
    Gap_S[j] = r_det[j] * dur_s[j];
  }
}
