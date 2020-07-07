data {
  int<lower=0> N; 
  int<lower=0> Asym;
  int<lower=0> Sym;
  int<lower=0> CS;
  
  int<lower=0> n_t;
  int<lower=0> Pop[n_t]; 
  int<lower=0> Noti[n_t]; 
  
  real<lower=0> Years[n_t];
  int year_survey;
  real<lower=0> r_sc_l;
  real<lower=0> r_sc_u;
  real<lower=0> r_death;
}
parameters {
  real<lower=0> del;
  real<lower=0, upper=1> prv0;
  real<lower=-0.15, upper=0.15> adr;
  real<lower=0> r_sym;
  real<lower=0> r_aware;
  real<lower=r_sc_l, upper = r_sc_u> r_sc;
}
transformed parameters {
  real<lower=0> ra = r_sc;
  real<lower=0> rs = r_sc + r_death;
  real<lower=0> rc = r_sc + r_death + 1 / del;
  
  real<lower=0> a0 = (rs + r_aware - adr) / r_sym;
  real<lower=0> c0 = r_aware / (rc - adr);
  real<lower=0, upper=1> pr_a;
  real<lower=0, upper=1> pr_s;
  real<lower=0, upper=1> pr_c;

  real<lower=0, upper=1> prv_a;
  real<lower=0, upper=1> prv_s;
  real<lower=0, upper=1> prv_c;

  
  vector<lower=0>[n_t] prv;
  vector<lower=0>[n_t] nr;
  

  pr_a = a0 / (a0 + 1 + c0);
  pr_s = 1 / (a0 + 1 + c0);
  pr_c = c0 / (a0 + 1 + c0);

  prv_a = prv0 * pr_a;
  prv_s = prv0 * pr_s;
  prv_c = prv0 * pr_c;

  
  for (i in 1:n_t) {
    prv[i] = prv0 * exp(- adr * (Years[i] - year_survey));
  }
  
  nr = prv * pr_c / del;
}
model {
  prv0 ~ uniform(0, 1);
  del ~ gamma(1E-1, 1E-1);

  r_aware ~ uniform(0, 24);
  r_sym ~ uniform(0, 24);
  r_sc ~ uniform(r_sc_l, r_sc_u);
  
  adr ~ uniform(-0.15, 0.15);

  target += binomial_lpmf(Asym | N, prv_a);
  target += binomial_lpmf(Sym | N, prv_s);
  
  for (i in 1:n_t) {
    target += poisson_lpmf(Noti[i] | nr[i] * Pop[i]);
  }
}
generated quantities {
  vector<lower=0>[n_t] inc_a;
  vector<lower=0>[n_t] inc_s;
  vector<lower=0>[n_t] inc_c;
  vector<lower=0>[n_t] noti;
  
  inc_a = (ra * pr_a + rs * pr_s + rc * pr_c - adr) * prv;
  inc_s = r_sym * pr_a * prv;
  inc_c = r_aware * pr_s * prv;

  noti = prv * pr_c / del;
}
