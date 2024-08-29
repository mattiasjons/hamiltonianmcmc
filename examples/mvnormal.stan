data {
  int<lower=0> N;               // number of observations
  int<lower=0> p;               // number of features
  matrix[N, p] X;               // predictor matrix
  vector[N] y;                  // outcome vector
}

parameters {
  vector[p] beta;               // coefficients for predictors
  //real<lower=0> sigma;          // error scale
}

model {
  beta ~ student_t(3, 0, 1);
  //beta ~ normal(0, 10);         // uninformative prior for coefficients
  //sigma ~ inv_gamma(3, 1);       // uninformative prior for variance
  y ~ normal(X * beta, 0.3696445);  // likelihood
}
