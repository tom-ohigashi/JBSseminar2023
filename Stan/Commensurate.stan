data {
  int<lower=0> x_h;
  int<lower=0> n_h;
  int<lower=0> x_CC;
  int<lower=0> n_CC;
  int<lower=0> x_CT;
  int<lower=0> n_CT;
  real<lower=0, upper=1> w;
  real<lower=0> lower_slab;
  real<lower=0> upper_slab;
  real<lower=0> spike;
}

parameters {
  real theta_CC;
  real theta_h;
  real<lower=0, upper=1> pi_CT;
  real<lower=0> tau;
}

transformed parameters{
  real<lower=0> inv_tau2;
  real<lower=0, upper=1> pi_CC = inv_logit(theta_CC);
  real<lower=0, upper=1> pi_h = inv_logit(theta_h);
  inv_tau2 = 1/(tau^2);
}

model {
  target += beta_lpdf(pi_h | 0.5, 0.5);
  target += log_sum_exp(log(w) + normal_lpdf(inv_tau2 | spike, 0.1), 
      log1m(w) + uniform_lpdf(inv_tau2 | lower_slab, upper_slab));
  target += normal_lpdf(theta_CC | theta_h, tau);
  target += binomial_lpmf(x_h | n_h, pi_h);
  target += binomial_lpmf(x_CC | n_CC, pi_CC);
  target += beta_lpdf(pi_CT | 0.5, 0.5);
  target += binomial_lpmf(x_CT | n_CT, pi_CT);
}

generated quantities {
  real g_pi = pi_CT - pi_CC;
}
