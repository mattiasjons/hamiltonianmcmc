data {
  int<lower=0> N; // number of observations
  vector[N] y;    // observations
}
parameters {
  real mu;        // mean parameter
  real<lower=0> sigma;  // standard deviation parameter
}
model {
  mu ~ normal(0, 3);
  sigma ~ inv_gamma(1, 1);
  y ~ normal(mu, 1);
}
