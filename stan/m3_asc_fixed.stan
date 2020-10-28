data {
  // Data from a prevalence survey
  int n_gp;
  int<lower=0> N[n_gp];
  int<lower=0> Asym[n_gp];
  int<lower=0> Sym[n_gp];
  int<lower=0> CS[n_gp];
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

  // Exogenous variables
  real<lower=0> r_death_bg[n_gp];
  real<lower=0> r_death_tb[n_gp];
}
parameters {
  real<lower=0> r_det;

  real<lower=0> r_sym;
  real<lower=0> r_aware;
  real<lower=r_sc_l, upper = r_sc_u> r_sc;
  
  vector<lower=0, upper=1>[n_gp] prv0;
  vector<lower=-0.15, upper=0.15>[n_gp] adr;  
}
transformed parameters {
  vector<lower=0>[n_gp] ra;
  vector<lower=0>[n_gp] rs;
  vector<lower=0>[n_gp] rc;

  vector<lower=0>[n_gp] a0;
  vector<lower=0>[n_gp] c0;
  
  vector<lower=0, upper=1>[n_gp] pr_a;
  vector<lower=0, upper=1>[n_gp] pr_s;
  vector<lower=0, upper=1>[n_gp] pr_c;
  vector<lower=0, upper=1>[n_gp] prv_a;
  vector<lower=0, upper=1>[n_gp] prv_s;
  vector<lower=0, upper=1>[n_gp] prv_c;

  matrix<lower=0>[n_gp, n_t] prv;
  matrix<lower=0>[n_gp, n_t] nr;
  
  
  for (j in 1:n_gp) {
    ra[j] = r_sc + r_death_bg[j];
    rs[j] = r_sc + r_death_bg[j] + r_death_tb[j];
    rc[j] = r_sc + r_death_bg[j] + r_death_tb[j];
    
    a0[j] = (rs[j] + r_aware - adr[j]) / r_sym;
    c0[j] = r_aware / (rc[j] + r_det - adr[j]);
    
    pr_a[j] = a0[j] / (a0[j] + 1 + c0[j]);
    pr_s[j] = 1 / (a0[j] + 1 + c0[j]);
    pr_c[j] = c0[j] / (a0[j] + 1 + c0[j]);
  
    prv_a[j] = prv0[j] * pr_a[j];
    prv_s[j] = prv0[j] * pr_s[j];
    prv_c[j] = prv0[j] * pr_c[j];
  }
  
  for (i in 1:n_t) {
    for (j in 1:n_gp) {
      prv[j, i] = prv0[j] * exp(- adr[j] * (Years[i] - YearSurveyed));
      nr[j, i] = prv[j, i] * pr_c[j] * r_det;
    }
  }
}
model {
  prv0 ~ uniform(0, 1);
  
  r_sym ~ inv_gamma(scale_dur, scale_dur);
  r_aware ~ inv_gamma(scale_dur, scale_dur);
  r_det ~ inv_gamma(scale_dur, scale_dur);
  
  r_sc ~ uniform(r_sc_l, r_sc_u);
  
  adr ~ uniform(-0.15, 0.15);

  for (j in 1:n_gp) {
    target += binomial_lpmf(Asym[j] | N[j], prv_a[j]);
    target += binomial_lpmf(Sym[j] | N[j], prv_s[j]);
    target += binomial_lpmf(CS[j] | N[j], prv_c[j]);
    
    for (i in 1:n_t) {
      target += poisson_lpmf(Noti[j, i] | nr[j, i] * Pop[j, i]);
    }
  }
}
generated quantities {
  matrix<lower=0>[n_gp, n_t] inc_a;
  matrix<lower=0>[n_gp, n_t] inc_s;
  matrix<lower=0>[n_gp, n_t] inc_c;
  matrix<lower=0>[n_gp, n_t] noti;
  
  real log_lik_noti[n_t, n_gp];
  real log_lik_pr[2, n_gp];
  
  real dur_a[n_gp];
  real dur_s[n_gp];
  real dur_c[n_gp];


  for (j in 1:n_gp) {
    // prevalence to prevalence survey data
    log_lik_pr[1, j] = binomial_lpmf(Asym[j] | N[j], prv_a[j]);
    log_lik_pr[2, j] = binomial_lpmf(Sym[j] | N[j], prv_s[j]);

    // notification rate to notification data
    for (i in 1:n_t) {
      log_lik_noti[i, j] = poisson_lpmf(Noti[j, i] | nr[j, i] * Pop[j, i]);
      
      inc_a[j, i] = (ra[j] * pr_a[j] + rs[j] * pr_s[j] + (rc[j] + r_det) * pr_c[j] - adr[j]) * prv[j, i];
      inc_s[j, i] = r_sym * pr_a[j] * prv[j, i];
      inc_c[j, i] = r_aware * pr_s[j] * prv[j, i];
    
      noti[j, i] = prv[j, i] * pr_c[j] * r_det;
    }
    
    dur_a[j] = 1 / (ra[j] + r_sym);
    dur_s[j] = 1 / (rs[j] + r_aware);
    dur_c[j] = 1 / (rc[j] + r_det);
  }
}
