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
  real<lower=0, upper=1> pr_fp;

  // Exogenous variables
  real<lower=0> r_death_a[n_gp];
  real<lower=0> r_death_sn[n_gp];
  real<lower=0> r_death_sp[n_gp];
}
parameters {
  real<lower=0, upper=1> p_sp; // Probability of sm+ at symptom onset
  real<lower=0> r_tr; // conversion rate, sm- to sm+
  real<lower=r_sc_l, upper = r_sc_u> r_sc; // self-cure rate
  
  vector<lower=0>[n_gp] r_det_sn; // delay to notification, sm-
  vector<lower=0>[n_gp] r_det_sp; // delay to notification, sm+
  vector<lower=0>[n_gp] r_sym; // symptom development rate == inversion of asymptomatic phase
  
  vector<lower=0, upper=1>[n_gp] prv0; // prevalence of all active TB
  vector<lower=-0.2, upper=0.2>[n_gp] adr; // annual decline rate
}
transformed parameters {
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

  

  for (j in 1:n_gp) {
    ra[j] = r_sc + r_death_a[j];
    rn[j] = r_sc + r_death_sn[j];
    rp[j] = r_sc + r_death_sp[j];
    
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
  r_tr ~ uniform(0, 2);
  r_sc ~ uniform(r_sc_l, r_sc_u);

  prv0 ~ uniform(0, 1);
  r_det_sp ~ inv_gamma(scale_dur, scale_dur);
  r_det_sn ~ inv_gamma(scale_dur, scale_dur);
  r_sym ~ inv_gamma(scale_dur, scale_dur);
  adr ~ uniform(-0.2, 0.2);
  
  for (j in 1:n_gp) {
    // prevalence to prevalence survey data
    target += binomial_lpmf(Asym[j] | N[j], prv_a[j]);
    target += binomial_lpmf(Sn[j] | N[j], prv_sn[j]);
    target += binomial_lpmf(Sp[j] | N[j], prv_sp[j]);

    // notification rate to notification data
    for (i in 1:n_t) {
      target += poisson_lpmf(NotiSn[j, i] | nr_sn[j, i] * Pop[j, i] / (1 - pr_fp));
      target += poisson_lpmf(NotiSp[j, i] | nr_sp[j, i] * Pop[j, i] / (1 - pr_fp));
    }
  }
}
generated quantities {
  real log_lik_noti[n_t, 2, n_gp];
  real log_lik_pr[3, n_gp];

  real dur_a[n_gp];
  real dur_sn[n_gp];
  real dur_sp[n_gp];
  
  real<lower=0> IncN_A[n_gp];
  real<lower=0> IncN_S[n_gp];
  real<lower=0> NotiN_Sn[n_gp];
  real<lower=0> NotiN_Sp[n_gp];
  real<lower=0> PrvN_A[n_gp];
  real<lower=0> PrvN_Sn[n_gp];
  real<lower=0> PrvN_Sp[n_gp];
  real<lower=0> CDR_A[n_gp];
  real<lower=0> CDR_S[n_gp];
  real<lower=0> CDR_Sp[n_gp];
  real<lower=0> Gap_A[n_gp];
  real<lower=0> Gap_S[n_gp];
  real<lower=0> Gap_Sp[n_gp];
  real<lower=0> Delay_Sn[n_gp];
  real<lower=0> Delay_Sp[n_gp];
  

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

    
    NotiN_Sn[j] = nr_sn[j, n_t] * Pop[j, n_t];
    NotiN_Sp[j] = nr_sp[j, n_t] * Pop[j, n_t];
    
    PrvN_A[j] = pr_a[j] * prv[j, n_t] * Pop[j, n_t];
    PrvN_Sn[j] = pr_sn[j] * prv[j, n_t] * Pop[j, n_t];
    PrvN_Sp[j] = pr_sp[j] * prv[j, n_t] * Pop[j, n_t];
    
    IncN_A[j] = (ra[j] * pr_a[j] + (rp[j] + r_det_sp[j]) * pr_sp[j] + (rn[j] + r_det_sn[j]) * pr_sn[j] - adr[j]) * prv[j, n_t] * Pop[j, n_t];
    IncN_S[j] = r_sym[j] * pr_a[j] * prv[j, n_t] * Pop[j, n_t];
    
    CDR_A[j] = (NotiN_Sp[j] + NotiN_Sn[j]) / IncN_A[j];
    CDR_S[j] = (NotiN_Sp[j] + NotiN_Sn[j]) / IncN_S[j];
    CDR_Sp[j] = NotiN_Sp[j] / (IncN_S[j] * p_sp);
    
    dur_a[j] = 1 / (ra[j] + r_sym[j]);
    dur_sn[j] = 1 / (rn[j] + r_tr + r_det_sn[j]);
    dur_sp[j] = 1 / (rp[j] + r_det_sp[j]);
    Delay_Sp[j] = dur_sp[j];
    Delay_Sn[j] = dur_sn[j] + r_tr * dur_sn[j]  * dur_sp[j];
    
    Gap_A[j] = r_sym[j] * dur_a[j];
    Gap_Sp[j] = r_det_sp[j] * dur_sp[j];
    Gap_S[j] = ((1 - p_sp) * (r_det_sn[j] + r_tr) * dur_sn[j] + p_sp) * Gap_Sp[j];
    
  }
}
