data {
  int<lower=0> N; 
  int<lower=0> Sn; 
  int<lower=0> Sp; 
  
  int<lower=0> n_t;
  int<lower=0> Pop[n_t]; 
  int<lower=0> Noti_Sp[n_t]; 
  int<lower=0> Noti_Sn[n_t]; 
  
  int<lower=0> Years[n_t];
  int year_survey;
  
  real<lower=0> r_death_sn;
  real<lower=0> r_death_sp;
  real<lower=0> r_sc_l;
  real<lower=0> r_sc_u;
}
parameters {
  real<lower=0, upper=1> p_sp;
  real<lower=0> del_sn;
  real<lower=0> del_sp;
  real<lower=0> r_tr;
  real<lower=0, upper=1> prv0;
  real<lower=-0.15, upper=0.15> adr;
  real<lower=r_sc_l, upper = r_sc_u> r_sc;
}
transformed parameters {
  real<lower=0> rn;
  real<lower=0> rp;
  
  real<lower=0> sn0;
  real<lower=0> sp0;
  real<lower=0, upper=1> pr_sn;
  real<lower=0, upper=1> pr_sp;
  real<lower=0, upper=1> prv_sn;
  real<lower=0, upper=1> prv_sp;
  
  vector<lower=0>[n_t] prv;
  vector<lower=0>[n_t] nr_sn;
  vector<lower=0>[n_t] nr_sp;
  vector<lower=0>[n_t] nr;
  
  rn = r_sc + r_death_sn + 1 / del_sn - adr;
  rp = r_sc + r_death_sp + 1 / del_sp - adr;
  
  sn0 = (1 - p_sp) * rp;    
  sp0 = p_sp * rn + r_tr;

  pr_sn = sn0 / (sn0 + sp0);
  pr_sp = sp0 / (sn0 + sp0);
  
  prv_sn = prv0 * pr_sn;
  prv_sp = prv0 * pr_sp;
  
  for (i in 1:n_t) {
    prv[i] = prv0 * exp(- adr * (Years[i] - year_survey));
  }
  
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
  
  r_sc ~ uniform(r_sc_l, r_sc_u);
  adr ~ uniform(-0.15, 0.15);

  target += binomial_lpmf(Sn | N, prv_sn);
  target += binomial_lpmf(Sp | N, prv_sp);
  
  for (i in 1:n_t) {
    target += poisson_lpmf(Noti_Sn[i] | nr_sn[i] * Pop[i]);
    target += poisson_lpmf(Noti_Sp[i] | nr_sp[i] * Pop[i]);
  }
}
generated quantities {
  vector<lower=0>[n_t] inc_sn;
  vector<lower=0>[n_t] inc_sp;
  vector<lower=0>[n_t] inc_s;
  vector<lower=0>[n_t] noti_sn;
  vector<lower=0>[n_t] noti_sp;
  
  inc_sn = (1 - p_sp) * (rn * pr_sn + rp * pr_sp - adr) * prv;
  inc_sp = (p_sp) * (rn * pr_sn + rp * pr_sp - adr) * prv;
  inc_s = inc_sn + inc_sp;
  noti_sn = prv * pr_sn / del_sn;
  noti_sp = prv * pr_sp / del_sp;
}