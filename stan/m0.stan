data {
  // Data from a prevalence survey
  int<lower=0> N;
  int<lower=0> Asym;
  int<lower=0> Sn;
  int<lower=0> Sp;
  real YearSurveyed; // timing of the survey

  // Time-series of notification data
  int<lower=0> n_t; // number of time points
  int<lower=0> Pop[n_t]; // population size
  int<lower=0> NotiSn[n_t]; // sm- counts
  int<lower=0> NotiSp[n_t]; // sm+ counts
  real<lower=0> Years[n_t]; // years of the notification data

  // Prior knowledge
  real<lower=0> r_sc_l;
  real<lower=0> r_sc_u;
  real<lower=0, upper=1> pr_fp;
  
  real<lower=0> r_death_a;
  real<lower=0> r_death_sn;
  real<lower=0> r_death_sp;
  real<lower=0> scale_dur;
}
parameters {
  real<lower=0, upper=1> p_sp; // Probability of sm+ at symptom onset
  real<lower=0> r_det_sn; // rate to notification, sm-
  real<lower=0> r_det_sp; // rate to notification, sm+
  real<lower=0> r_tr; // conversion rate, sm- to sm+
  real<lower=0, upper=1> prv0; // prevalence of all active TB
  real<lower=-0.2, upper=0.2> adr; // annual decline rate
  real<lower=0> r_sym; // symptom development rate == inversion of asymptomatic phase
  real<lower=r_sc_l, upper = r_sc_u> r_sc; // self-cure rate
}
transformed parameters {
  real<lower=0> ra = r_sc + r_death_a;
  real<lower=0> rn = r_sc + r_death_sn;
  real<lower=0> rp = r_sc + r_death_sp;

  real<lower=0> sn0;
  real<lower=0> sp0;
  real<lower=0, upper=1> pr_a;
  real<lower=0, upper=1> pr_sn;
  real<lower=0, upper=1> pr_sp;
  real<lower=0, upper=1> prv_a;
  real<lower=0, upper=1> prv_sn;
  real<lower=0, upper=1> prv_sp;

  vector<lower=0>[n_t] prv;
  vector<lower=0>[n_t] nr_sn;
  vector<lower=0>[n_t] nr_sp;
  vector<lower=0>[n_t] nr;


  sn0 = (1 - p_sp) * r_sym / (r_tr + rn + r_det_sn - adr);
  sp0 = (p_sp * r_sym + r_tr * sn0) / (rp + r_det_sp - adr);

  // probabilities of the states among all active TB
  pr_a = 1 / (1 + sn0 + sp0);
  pr_sn = sn0 * pr_a;
  pr_sp = sp0 * pr_a;

  // prevalence of the states
  prv_a = prv0 * pr_a;
  prv_sn = prv0 * pr_sn;
  prv_sp = prv0 * pr_sp;


  // forecasts/backcasts of prevalence
  for (i in 1:n_t) {
    prv[i] = prv0 * exp(- adr * (Years[i] - YearSurveyed));
  }

  // notification rates
  nr_sn = prv * pr_sn * r_det_sn;
  nr_sp = prv * pr_sp * r_det_sp;
  nr = nr_sn + nr_sp;
}
model {
  p_sp ~ uniform(0, 1);
  r_tr ~ uniform(0, 2);
  r_sc ~ uniform(r_sc_l, r_sc_u);

  r_det_sp ~ inv_gamma(scale_dur, scale_dur);
  r_det_sn ~ inv_gamma(scale_dur, scale_dur);
  r_sym ~ inv_gamma(scale_dur, scale_dur);
  adr ~ uniform(-0.2, 0.2);

  prv0 ~ uniform(0, 1);


  // prevalence to prevalence survey data
  target += binomial_lpmf(Asym | N, prv_a);
  target += binomial_lpmf(Sn | N, prv_sn);
  target += binomial_lpmf(Sp | N, prv_sp);

  // notification rate to notification data
  for (i in 1:n_t) {
    target += poisson_lpmf(NotiSn[i] | nr_sn[i] * Pop[i] / (1 - pr_fp));
    target += poisson_lpmf(NotiSp[i] | nr_sp[i] * Pop[i] / (1 - pr_fp));
  }
}
generated quantities {
  real dur_a;
  real dur_sn;
  real dur_sp;

  
  real<lower=0> IncN_A;
  real<lower=0> IncN_S;
  real<lower=0> NotiN_Sn;
  real<lower=0> NotiN_Sp;
  real<lower=0> PrvN_A;
  real<lower=0> PrvN_Sn;
  real<lower=0> PrvN_Sp;
  real<lower=0> CDR_A;
  real<lower=0> CDR_S;
  real<lower=0> CDR_Sp;
  real<lower=0> Gap_A;
  real<lower=0> Gap_S;
  real<lower=0> Gap_Sp;
  real<lower=0> Delay_Sn;
  real<lower=0> Delay_Sp;
  
  
  dur_a = 1 / (ra + r_sym);
  dur_sn = 1 / (rn + r_tr + r_det_sn);
  dur_sp = 1 / (rp + r_det_sp);

  NotiN_Sn = nr_sn[n_t] * Pop[n_t];
  NotiN_Sp = nr_sp[n_t] * Pop[n_t];
  
  PrvN_A = pr_a * prv[n_t] * Pop[n_t];
  PrvN_Sn = pr_sn * prv[n_t] * Pop[n_t];
  PrvN_Sp = pr_sp * prv[n_t] * Pop[n_t];
  
  IncN_A = (ra * pr_a + (rp + r_det_sp) * pr_sp + (rn + r_det_sn) * pr_sn - adr) * prv[n_t] * Pop[n_t];
  IncN_S = r_sym * pr_a * prv[n_t] * Pop[n_t];
  
  CDR_A = (NotiN_Sp + NotiN_Sn) / IncN_A;
  CDR_S = (NotiN_Sp + NotiN_Sn) / IncN_S;
  CDR_Sp = NotiN_Sp / (IncN_S * p_sp);
  
  dur_a = 1 / (ra + r_sym);
  dur_sn = 1 / (rn + r_tr + r_det_sn);
  dur_sp = 1 / (rp + r_det_sp);
  Delay_Sp = dur_sp;
  Delay_Sn = dur_sn + r_tr * dur_sn  * dur_sp;
  
  Gap_A = r_sym * dur_a;
  Gap_Sp = r_det_sp * dur_sp;
  Gap_S = ((1 - p_sp) * (r_det_sn + r_tr) * dur_sn + p_sp) * Gap_Sp;

}
