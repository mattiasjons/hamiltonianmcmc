# hamiltonianmcmc

R package to run Hamiltonian MCMC sampling

# How to Install

Easiest way to install is:

```r
devtools::install_github("mattiasjons/hamiltonianmcmc", upgrade_dependencies = FALSE)
```

# How to use
```r
# Function to calculate the negative log likelihood (up to a constant)
negative_log_likelihood <- function(x, mean, covariance) {
  k <- length(mean)
  det_cov <- det(covariance)
  stopifnot("Covariance matrix must be positive definite."= det_cov > 0)

  constant <- -0.5 * (k * log(2 * pi) + log(det_cov))
  exponent <- 0.5 * t(x - mean) %*% solve(covariance) %*% (x - mean)
  log_likelihood <- constant + exponent

  return(log_likelihood)
}

# Function to calculate the gradient of the negative log likelihood
gradient_negative_log_likelihood <- function(x, mean, covariance) {
  k <- length(mean)
  det_cov <- det(covariance)
  stopifnot("Covariance matrix must be positive definite."= det_cov > 0)

  gradient <- -solve(covariance) %*% (x - mean)

  return(gradient)
}

set.seed(123)
mean_vector <- c(2, -2)
covariance_matrix <- matrix(c(1, 0.98, 0.98, 2), nrow = 2, ncol = 2)

num_samples <- 1000
initial_state <- c(0, 0)
step_size <- 0.05
num_steps <- 30

# Define the multivariate normal density and gradient functions
density_function <- function(x) -negative_log_likelihood(x, mean_vector, covariance_matrix)
gradient_function <- function(x) -gradient_negative_log_likelihood(x, mean_vector, covariance_matrix)

samples <- hamiltonian_mcmc(initial_state, num_samples,
                            step_size, num_steps,
                            density_function, gradient_function,
                            M=c(1, 1))
plot(samples[,1])
```
