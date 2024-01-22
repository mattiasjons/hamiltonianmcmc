#' Hamiltonian Markov Chain Monte Carlo sampler.
#'
#' @description
#'
#' `hamiltonian_mcmc()` is the main function of this library, that provides a simple way to run a Hamiltonian MCMC sampler for given density and gradient functions.
#'
#' At the moment, no vignette exists.
#

#' Leapfrog integrator
#'
#' @param current_state A number.
#' @param momentum A number.
#' @param step_size A number.
#' @param gradient_function A number.
#' @param num_steps A number.
#' @export
#' @returns A list of proposed state and proposed momentum
leapfrog_integration <- function(current_state, momentum, step_size, gradient_function, num_steps) {
  proposed_state <- current_state

  proposed_momentum <- momentum - 0.5 * step_size * gradient_function(current_state)

  for (j in 1:num_steps) {
    proposed_state <- proposed_state + step_size * proposed_momentum
    proposed_momentum <- proposed_momentum - step_size * gradient_function(proposed_state)
  }

  proposed_momentum <- proposed_momentum - 0.5 * step_size * gradient_function(proposed_state)

  return(list(proposed_state = proposed_state, proposed_momentum = proposed_momentum))
}


#' Sample using hamiltonian MCMC
#'
#' @param initial_state A vector with length \eqn{d}, representing the initial state.
#' @param num_samples Number of samples to generate.
#' @param step_size How large is the step size in the leapfrof integration.
#' @param num_steps Number of leap frog steps. If Adaptive Mode is enabled, this represents the initial value.
#' @param density_function A function representing the target density to sample from. The function should take a single vector valued parameter of length \eqn{d}.
#' @param gradient_function A function representing the gradient of the density function. The function should take a single vector valued parameter of length \eqn{d}.
#' @param M The mass matrix for sampling momenta. Can be passed as either a full \eqn{d \times d} matrix or as a vector of length \eqn{d}, interpreted as the diagnoal of the mass matrix
#' @param adaptive_mode Boolean to indicate whether the step size should be dynamic based on a target acceptance rate. FALSE by default.
#' @param adaptive_target_acceptance If adaptive mode is enabled, this sets the target acceptance rate. Default = \eqn{0.85}
#' @returns Returns a matrix of samples.
#' @export
#' @examples
#'# Function to calculate the negative log likelihood (up to a constant)
#'negative_log_likelihood <- function(x, mean, covariance) {
#'  k <- length(mean)
#'  det_cov <- det(covariance)
#'  stopifnot("Covariance matrix must be positive definite."= det_cov > 0)
#'
#'  constant <- -0.5 * (k * log(2 * pi) + log(det_cov))
#'  exponent <- 0.5 * t(x - mean) %*% solve(covariance) %*% (x - mean)
#'  log_likelihood <- constant + exponent
#'
#'  return(log_likelihood)
#'}
#'
#'# Function to calculate the gradient of the negative log likelihood
#'gradient_negative_log_likelihood <- function(x, mean, covariance) {
#'  k <- length(mean)
#'  det_cov <- det(covariance)
#'  stopifnot("Covariance matrix must be positive definite."= det_cov > 0)
#'
#'  gradient <- -solve(covariance) %*% (x - mean)
#'
#'  return(gradient)
#'}
#'
#'set.seed(123)
#'mean_vector <- c(2, -2)
#'covariance_matrix <- matrix(c(1, 0.98, 0.98, 2), nrow = 2, ncol = 2)

#'num_samples <- 1000
#'initial_state <- c(0, 0)
#'step_size <- 0.05
#'num_steps <- 30

#'# Define the multivariate normal density and gradient functions
#'density_function <- function(x) -negative_log_likelihood(x, mean_vector, covariance_matrix)
#'gradient_function <- function(x) -gradient_negative_log_likelihood(x, mean_vector, covariance_matrix)

#'samples <- hamiltonian_mcmc(initial_state, num_samples,
#'                            step_size, num_steps,
#'                            density_function, gradient_function,
#'                            M=c(1, 1))
#'plot(samples[,1])
hamiltonian_mcmc <- function(initial_state,
                             num_samples,
                             step_size,
                             num_steps,
                             density_function,
                             gradient_function,
                             M,
                             adaptive_mode = FALSE,
                             adaptive_target_acceptance = 0.85) {
  #TODO: Asserts

  if (length(M)==length(initial_state) & length(M)>1) {
    M = diag(M)
  }

  current_state <- initial_state
  samples <- matrix(nrow = num_samples, ncol = length(initial_state))

  for (i in 1:num_samples) {
    momentum <- MASS::mvrnorm(1, mu = rep(0, length(current_state)), Sigma = M)

    # Leapfrog integration
    leapfrog_result <- leapfrog_integration(current_state, momentum, step_size, gradient_function, num_steps)
    proposed_state <- leapfrog_result$proposed_state
    proposed_momentum <- leapfrog_result$proposed_momentum

    # Metropolis acceptance step
    current_energy <- -density_function(current_state) + 0.5 * sum(momentum^2)
    proposed_energy <- -density_function(proposed_state) + 0.5 * sum(proposed_momentum^2)

    acceptance_ratio <- exp(current_energy - proposed_energy)

    if (stats::runif(1) < acceptance_ratio) {
      current_state <- proposed_state
      samples[i,] <- current_state
    }

    if (adaptive_mode) {
      if (acceptance_ratio > adaptive_target_acceptance) {
        step_size = step_size * 1.1
      } else {
        step_size = step_size / 2
      }
    }
  }

  return(samples)
}
