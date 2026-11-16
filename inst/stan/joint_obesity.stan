// Joint model: latent LAD-level obesity in early pregnancy + Poisson stillbirth model.
//
// Obesity sub-model (percent scale):
//   o_i = a0 + Z_i * gamma + u_{g(i)},  u_g ~ N(0, tau^2)      (g = upper-tier LA of LAD i)
// Observed at upper-tier level (births-weighted average of member LADs), with binomial
// sampling error and a period offset (2022/23 vs 2023/24):
//   zobs_k ~ N( sum_i W[k,i] o_i + delta * is_2223_k , se_k^2 + omega^2 )
// Optional LAD/provider-level observations (e.g. hospital booking data) enter the same way
// with W rows that pick out a single LAD.
//
// Stillbirth sub-model:
//   y_i ~ Poisson( B_i exp(beta0 + X_i beta + beta_ob * (o_i - ob_c)/ob_s + b_i) ),  b_i ~ N(0, sigma^2)
// Areas with in_lik == 0 are held out (cross-validation).
data {
  int<lower=1> N;                  // LADs
  int<lower=1> G;                  // upper-tier authorities
  array[N] int<lower=1, upper=G> grp;
  array[N] int<lower=0> y;
  vector[N] log_offset;
  int<lower=0> K;                  // stillbirth covariates (standardised), excluding obesity
  matrix[N, K] X;
  int<lower=1> P;                  // obesity-model covariates (standardised)
  matrix[N, P] Z;
  int<lower=0> M;                  // obesity observations
  matrix[M, N] Wobs;               // aggregation weights (rows sum to 1)
  vector[M] zobs;                  // observed % obese at booking
  vector<lower=0>[M] se_obs;       // binomial SE (percentage points)
  vector<lower=0, upper=1>[M] is_2223;
  real ob_c;                       // centring / scaling of obesity in the outcome model
  real<lower=0> ob_s;
  array[N] int<lower=0, upper=1> in_lik;
  real<lower=0> sigma_scale;      // half-normal prior scale for the random-effect SD
}
transformed data {
  int n_lik = sum(in_lik);
  array[n_lik] int idx;
  {
    int j = 1;
    for (i in 1:N) if (in_lik[i] == 1) { idx[j] = i; j += 1; }
  }
}
parameters {
  // obesity sub-model
  real a0;
  vector[P] gamma;
  real<lower=0> tau;
  vector[G] u_raw;
  real delta;
  real<lower=0> omega;
  // stillbirth sub-model
  real beta0;
  vector[K] beta;
  real beta_ob;
  real<lower=0> sigma;
  vector[N] theta;
}
transformed parameters {
  vector[N] ob = a0 + Z * gamma + tau * u_raw[grp];      // latent % obese in pregnancy
  vector[N] ob_std = (ob - ob_c) / ob_s;
  vector[N] b = sigma * theta;
}
model {
  // obesity sub-model
  a0 ~ normal(26, 10);
  gamma ~ normal(0, 5);
  tau ~ normal(0, 5);
  u_raw ~ std_normal();
  delta ~ normal(0, 2);
  omega ~ normal(0, 2);
  zobs ~ normal(Wobs * ob + delta * is_2223, sqrt(square(se_obs) + square(omega)));

  // stillbirth sub-model
  {
    vector[N] eta = log_offset + beta0 + X * beta + beta_ob * ob_std + b;
    y[idx] ~ poisson_log(eta[idx]);
  }
  beta0 ~ normal(-5.5, 2);
  beta ~ normal(0, 0.5);
  beta_ob ~ normal(0, 0.5);
  sigma ~ normal(0, sigma_scale);  // half-normal prior on the random-effect SD
  theta ~ std_normal();
}
generated quantities {
  vector[N] log_rate = beta0 + X * beta + beta_ob * ob_std + b;
  vector[N] log_lik;
  for (i in 1:N) log_lik[i] = poisson_log_lpmf(y[i] | log_offset[i] + log_rate[i]);
}
