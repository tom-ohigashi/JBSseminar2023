data {
  int<lower=0> H;
  array[H] int<lower=0> x_h;
  array[H] int<lower=0> n_h;
  int<lower=0> x_CC;
  int<lower=0> n_CC;
  int<lower=0> x_CT;
  int<lower=0> n_CT;
  real<lower=0> betascale;
  int<lower=0> nu;
}

parameters {
  real<lower=0> tau;
  vector<lower=0>[H] lambda;
  vector[H] beta_raw;
  real theta_CC;
  real<lower=0, upper=1> pi_CT;
}

transformed parameters{
  vector[H] theta_h;
  vector[H] beta;
  beta = tau * lambda .* beta_raw;
  theta_h = theta_CC + beta;
}

model {
  target += normal_lpdf(theta_CC | 0, 100);
  target += normal_lpdf(beta_raw | 0, 1);
  target += student_t_lpdf(tau | nu, 0, betascale);
  target += student_t_lpdf(lambda | 1, 0, 1);
  target += binomial_logit_lpmf(x_h | n_h, theta_h);
  target += binomial_logit_lpmf(x_CC | n_CC, theta_CC);
  target += beta_lpdf(pi_CT | 0.5, 0.5);
  target += binomial_lpmf(x_CT | n_CT, pi_CT);
}

generated quantities {
  real pi_CC = inv_logit(theta_CC);
  real g_pi = pi_CT - pi_CC;
}
