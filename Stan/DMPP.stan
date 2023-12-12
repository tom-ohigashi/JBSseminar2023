data {
  int<lower=0> H;
  array[H] int<lower=0> x_h;
  array[H] int<lower=0> n_h;
  int<lower=0> x_CC;
  int<lower=0> n_CC;
  int<lower=0> x_CT;
  int<lower=0> n_CT;
}

parameters {
  real<lower=0, upper=1> pi_CC;
  real<lower=0, upper=1> pi_CT;
  array[H] real<lower=0, upper=1> delta;
  real<lower=0> a;
  real<lower=0> b;
}

transformed parameters{
  real<lower=0, upper=1> mu;
  real<lower=0> sigmasq;
  real lik_CC;
  real SC;
  array[H] real delta_sum1;
  array[H] real delta_sum2;
  array[H] real delta_sum3;

  mu = a / (a + b);
  sigmasq = (mu * (1 - mu)) / (a + b + 1);

  for(h in 1:H){
    delta_sum1[h] = delta[h] * x_h[h];
    delta_sum2[h] = delta[h] * (n_h[h] - x_h[h]);
    delta_sum3[h] = delta[h] * n_h[h];
  }
  
  lik_CC = (sum(delta_sum1) + x_CC + 0.5 - 1) * log(pi_CC) + 
          (sum(delta_sum2) + (n_CC - x_CC) + 0.5 - 1) * log(1 - pi_CC);
  SC = lgamma(sum(delta_sum1) + 0.5) + lgamma(sum(delta_sum2) + 0.5) - 
        lgamma(sum(delta_sum3) + 0.5 + 0.5);
}

model {
  target += uniform_lpdf(mu | 0, 1);
  target += inv_gamma_lpdf(sigmasq | 1, 1);
  target += beta_lpdf(delta | a, b);
  target += lik_CC;
  target += -SC;
  target += beta_lpdf(pi_CT | 0.5, 0.5);
  target += binomial_lpmf(x_CT | n_CT, pi_CT);
}

generated quantities {
  real g_pi = pi_CT - pi_CC;
}
