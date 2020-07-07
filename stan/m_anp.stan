data {
  // Data from a prevalence survey
  int<lower=0> N; 
  int<lower=0> Asym;
  int<lower=0> Sn; 
  int<lower=0> Sp; 
  int year_survey; // timing of the survey
  
  // Time-series of notification data
  int<lower=0> n_t; // number of time points
  int<lower=0> Pop[n_t]; // population size
  int<lower=0> Noti_Sn[n_t]; // sm- counts
  int<lower=0> Noti_Sp[n_t]; // sm+ counts
  real<lower=0> Years[n_t]; // years of the notification data
  
  // Prior knowledge
  real<lower=0> r_sc_l;
  real<lower=0> r_sc_u;
  real<lower=0> r_death_sn;
  real<lower=0> r_death_sp;
}
parameters {
  real<lower=0, upper=1> p_sp; // Probability of sm+ at symptom onset
  real<lower=0> del_sn; // delay to notification, sm-
  real<lower=0> del_sp; // delay to notification, sm+
  real<lower=0> r_tr; // conversion rate, sm- to sm+
  real<lower=0, upper=1> prv0; // prevalence of all active TB
  real<lower=-0.15, upper=0.15> adr; // annual decline rate
  real<lower=0> r_sym; // symptom development rate == inversion of asymptomatic phase
  real<lower=r_sc_l, upper = r_sc_u> r_sc; // self-cure rate
}
transformed parameters {
  real<lower=0> ra = r_sc;
  real<lower=0> rn = r_sc + r_death_sn + 1 / del_sn;
  real<lower=0> rp = r_sc + r_death_sp + 1 / del_sp;
  
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
  

  sn0 = (1 - p_sp) * r_sym / (r_tr + rn - adr);    
  sp0 = (p_sp * r_sym + r_tr * sn0) / (rp - adr);
    
  // probabilities of the states among all active TB
  pr_a = 1 / (1 + sn0 + sp0);
  pr_sn = sn0 * pr_a;
  pr_sp = sp0 * pr_a;
  
  // prevalences of the states
  prv_a = prv0 * pr_a;
  prv_sn = prv0 * pr_sn;
  prv_sp = prv0 * pr_sp;
  
  
  // forecasts/backcasts of prevalences
  for (i in 1:n_t) {
    prv[i] = prv0 * exp(- adr * (Years[i] - year_survey));
  }
  
  // notification rates
  nr_sn = prv * pr_sn / del_sn;
  nr_sp = prv * pr_sp / del_sp;
  nr = nr_sn + nr_sp;
}
model {
  prv0 ~ uniform(0, 1);
  p_sp ~ uniform(0, 1);
  r_tr ~ uniform(0, 0.5);
  del_sp ~ gamma(1E-1, 1E-1);
  del_sn ~ gamma(1E-1, 1E-1);
  
  r_sym ~ uniform(0, 24);
  r_sc ~ uniform(r_sc_l, r_sc_u);
  
  adr ~ uniform(-0.15, 0.15);

  // prevalence to prevalence survey data
  target += binomial_lpmf(Asym | N, prv_a);
  target += binomial_lpmf(Sn | N, prv_sn);
  target += binomial_lpmf(Sp | N, prv_sp);
  
  // notification rate to notification data
  for (i in 1:n_t) {
    target += poisson_lpmf(Noti_Sn[i] | nr_sn[i] * Pop[i]);
    target += poisson_lpmf(Noti_Sp[i] | nr_sp[i] * Pop[i]);
  }
}
generated quantities {
  vector<lower=0>[n_t] inc_a;
  vector<lower=0>[n_t] inc_sn;
  vector<lower=0>[n_t] inc_sp;
  vector<lower=0>[n_t] inc_s;
  vector<lower=0>[n_t] noti_sn;
  vector<lower=0>[n_t] noti_sp;
  vector<lower=0>[n_t] ni;
  
  // incidence estimates
  inc_a = (ra * pr_a + rn * pr_sn + rp * pr_sp - adr) * prv;
  inc_sn = (1 - p_sp) * r_sym * pr_a * prv;
  inc_sp = (p_sp) * r_sym * pr_a * prv + r_tr * pr_sn * prv;
  inc_s = r_sym * pr_a * prv;
  noti_sn = prv * pr_sn / del_sn;
  noti_sp = prv * pr_sp / del_sp;

  for (i in 1:n_t) {
    ni[i] = (noti_sn[i] + noti_sp[i]) / inc_a[i];
  }
}
