// Poisson areal model for stillbirth counts with three options for the
// area-level random effect:
//   re_type = 0 : IID (Poisson-lognormal, no spatial structure)
//   re_type = 1 : BYM2 (Riebler et al. 2016), scaled ICAR + IID, non-centred
//   re_type = 2 : Leroux CAR (as used in CARBayes::S.CARleroux)
// Areas with in_lik == 0 are held out (for K-fold CV): their counts do not enter
// the likelihood but their random effects and predictions are still generated.
functions {
  real icar_lpdf(vector phi, array[] int node1, array[] int node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]);
  }
}
data {
  int<lower=1> N;
  int<lower=0> N_edges;
  array[N_edges] int<lower=1, upper=N> node1;
  array[N_edges] int<lower=1, upper=N> node2;
  array[N] int<lower=0> y;
  vector[N] log_offset;           // log(total births)
  int<lower=0> K;
  matrix[N, K] X;                 // standardised covariates
  array[N] int<lower=0, upper=1> in_lik;
  int<lower=0, upper=2> re_type;
  real<lower=0> scaling_factor;   // BYM2 ICAR scaling
  vector[N] lambda_W;             // eigenvalues of D - W (Leroux)
  matrix[re_type == 2 ? N : 0, re_type == 2 ? N : 0] V_W;   // eigenvectors of D - W (Leroux)
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
  real beta0;
  vector[K] beta;
  real<lower=0> sigma;
  real<lower=0, upper=1> rho;
  vector[N] theta;                // standard-normal latent (IID part / Leroux whitened field)
  vector[re_type == 1 ? N : 0] phi;   // ICAR component (BYM2 only)
}
transformed parameters {
  vector[N] b;                    // area random effect on the log-rate scale
  if (re_type == 0) {
    b = sigma * theta;
  } else if (re_type == 1) {
    b = sigma * (sqrt(rho / scaling_factor) * phi + sqrt(1 - rho) * theta);
  } else {
    // Leroux, non-centred: Q(rho) = rho (D - W) + (1 - rho) I = V diag(rho*lambda + 1 - rho) V'
    b = sigma * (V_W * (theta ./ sqrt(rho * lambda_W + 1 - rho)));
  }
}
model {
  vector[N] eta = log_offset + beta0 + X * beta + b;
  y[idx] ~ poisson_log(eta[idx]);

  beta0 ~ normal(-5.5, 2);
  beta ~ normal(0, 0.5);          // log-RR per SD: 95% prior within RR 0.38-2.7
  sigma ~ normal(0, sigma_scale);  // half-normal prior on the random-effect SD
  rho ~ beta(1, 1);

  theta ~ std_normal();
  if (re_type == 1) {
    phi ~ icar(node1, node2);
    sum(phi) ~ normal(0, 0.001 * N);
  }
}
generated quantities {
  vector[N] log_rate = beta0 + X * beta + b;     // log stillbirths per total birth
  vector[N] fixed_part = X * beta;
  vector[N] log_lik;
  array[N] int y_rep;
  for (i in 1:N) {
    real eta = log_offset[i] + log_rate[i];
    log_lik[i] = poisson_log_lpmf(y[i] | eta);
    y_rep[i] = poisson_log_rng(fmin(eta, 20));
  }
}
