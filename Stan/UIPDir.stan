data {
  int<lower=0> H;
  array[H] int<lower=0> x_h;
  array[H] int<lower=0> n_h;
  int<lower=0> x_CC;
  int<lower=0> n_CC;
  int<lower=0> x_CT;
  int<lower=0> n_CT;
  vector<lower=0>[H] gamma_h;
  vector<lower=0, upper=1>[H] theta_h;
  vector<lower=0>[H] I_U;
}

parameters {
  real<lower=0, upper=1> pi_CC;
  real<lower=0, upper=1> pi_CT;
  
  real<lower=0, upper=sum(n_h)> M;
  simplex[H] omega; 
}

transformed parameters{
  real<lower=0, upper=1> mu;
  real<lower=0> eta2;
  real<lower=0> alpha;
  real<lower=0> beta;

  mu = sum(omega .* theta_h);
  eta2 = 1/(M * sum(omega .* I_U));
  
  if(mu*(1-mu) > eta2){
    alpha = mu * ((mu*(1-mu)/eta2) - 1);
    beta = (1-mu) * ((mu*(1-mu)/eta2) - 1);
  }else{
    alpha = mu * 0.01;
    beta = (1-mu) * 0.01;
  }
}

model {
  target += dirichlet_lpdf(omega | gamma_h);
  target += uniform_lpdf(M | 0, sum(n_h));
  target += beta_lpdf(pi_CC | alpha, beta);
  target += binomial_lpmf(x_CC | n_CC, pi_CC);
  target += beta_lpdf(pi_CT | 0.5, 0.5);
  target += binomial_lpmf(x_CT | n_CT, pi_CT);
}

generated quantities {
  real g_pi = pi_CT - pi_CC;
  real theta_CC = logit(pi_CC);
}
