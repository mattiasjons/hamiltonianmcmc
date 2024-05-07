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
#' @returns A list of proposed state, proposed momentum and final gradient
leapfrog_integration <- function(current_state, momentum, step_size, tgt_density, num_steps) {

  proposed_state <- current_state
  proposed_momentum <- momentum - 0.5 * step_size * as.numeric(tgt_density$grad(current_state)[[1]])

  for (j in 1:num_steps) {
    proposed_state <- proposed_state - step_size * proposed_momentum
    proposed_momentum <- proposed_momentum - step_size * as.numeric(tgt_density$grad(proposed_state)[[1]])
  }

  final_grad <- as.numeric(tgt_density$grad(proposed_state)[[1]])
  proposed_momentum <- proposed_momentum - 0.5 * step_size * final_grad

  return(list(proposed_state = proposed_state, proposed_momentum = proposed_momentum,
              gradient = final_grad))
}


#' Sample using hamiltonian MCMC
#'
#' @param initial_state A vector with length \eqn{d}, representing the initial state.
#' @param num_samples Number of samples to generate.
#' @param step_size How large is the step size in the leapfrof integration.
#' @param num_steps Number of leap frog steps. If Adaptive Mode is enabled, this represents the initial value.
#' @param tgt_density An object of class TargetDensity representing the target density to sample from.
#' @param auto_mass_matrix A boolean indicating wether the Mass matrix should be adjusted dynamically based on accepted posterior samples
#' @param n_mass_matrix_comp If the mass matrix is dynamically adapted, it is done so using a eigendecomposition. This parameter sets the number of components for the decomposition.
#' @param M The mass matrix for sampling momenta. Can be passed as either a full \eqn{d \times d} matrix or as a vector of length \eqn{d}, interpreted as the diagnoal of the mass matrix
#' @param adaptive_mode Boolean to indicate whether the step size should be dynamic based on a target acceptance rate. FALSE by default.
#' @param adaptive_target_acceptance If adaptive mode is enabled, this sets the target acceptance rate. Default = \eqn{0.85}
#' @return A list of the following
#' \itemize{
#'   \item samples: A matrix of accepted samples. Each non accepted iteration will have NA entries.
#'   \item proposals: A matrix of all proposals, no matter if they were accepted or not. Used primarily for diagnostics.
#'   \item energy: A matrix with two columns. Column1 is the total hamiltonian for the state at each iteration, Column2 is the hamiltonian for the proposal generated in each iteration.
#'   \item metropolis_acceptance: A matrix with a single column, containing the metropolis acceptance (i.e. random $U(0,1)$) generated at each iteration.
#' }
#' @export
#' @examples
#' library(hamiltonianmcmc)
#' library(torch)
#' library(MASS)
#'
#' d = 3
#' rho <- runif(3)
#' sigma <- c(3, 2, 1)
#' df <- mvrnorm(100, rep(0, d), t(rho%*%diag(sigma))%*%t(rho))
#'
#' tgt_dens <- hamiltonianmcmc::TargetDensity$new(data = torch_tensor(df),
#'                                                init_mu = rnorm(d),
#'                                                diag(d))
#'
#' num_samples <- 1500
#' initial_state <- rnorm(d)
#' step_size <- 0.0005
#' num_steps <- 30
#'
#' hmc_res <- hamiltonian_mcmc(initial_state = initial_state, num_samples = num_samples,
#'                             step_size = step_size, num_steps = num_steps,
#'                             tgt_density = tgt_dens, auto_mass_matrix = F,
#'                             n_mass_matrix_comp = 1, M=rep(0.01, d),
#'                             adaptive_mode = T, adaptive_target_acceptance = 0.8)
#'
#' plot(hmc_res$samples[100:1500,])
hamiltonian_mcmc <- function(initial_state,
                             num_samples,
                             step_size,
                             num_steps,
                             tgt_density,
                             auto_mass_matrix,
                             n_mass_matrix_comp,
                             M,
                             adaptive_mode = FALSE,
                             adaptive_target_acceptance = 0.85) {
  #TODO: Asserts

  if(is.null(M)) {
    M <- rep(1, length(initial_state))
  }

  if (length(M)==length(initial_state) & length(M)>1) {
    M = diag(M)
  }

  if (auto_mass_matrix) {
    library(incrementalpca)
  }

  current_state <- initial_state
  samples <- matrix(nrow = num_samples, ncol = length(initial_state))
  proposals <- matrix(nrow=num_samples, ncol = length(initial_state))
  energy <- matrix(nrow = num_samples, ncol=2)
  metropolis_acceptance = matrix(nrow = num_samples, ncol=1)
  init_done <- F
  n_accepted <- 0
  last_updated_mass_step <- 0

  for (i in 1:num_samples) {
    if (i%%50==0) {
      cat(paste0('Generated ', i, ' samples.'))
      cat('\r\n')
    }

    if (auto_mass_matrix){
      if (!init_done && n_accepted > n_mass_matrix_comp) {
        M_decomp <- IncrementalDecomposition$new(samples[!is.na(samples[,1]),], n_mass_matrix_comp)
        init_done <- T
        last_updated_mass_step <- i
        M <- M_decomp$get_precision()
      } else if (init_done && i >= last_updated_mass_step + 5) {
        smp <- samples[(i-4):i,]
        M_decomp$partial_fit(smp, 0.98)
        last_updated_mass_step <- i
        M <- M_decomp$get_precision()
      }
    }

    momentum <- MASS::mvrnorm(1, mu = rep(0, length(current_state)), Sigma = M)

    # Leapfrog integration
    leapfrog_result <- leapfrog_integration(current_state, momentum, step_size, tgt_density, num_steps)
    proposed_state <- leapfrog_result$proposed_state
    proposed_momentum <- leapfrog_result$proposed_momentum

    # Metropolis acceptance step
    current_energy <- -as.numeric(tgt_density$value(current_state)) + 0.5 * t(momentum) %*% M %*% momentum
    energy[i,1] <- current_energy
    proposed_energy <- -as.numeric(tgt_density$value(proposed_state)) + 0.5 * t(momentum) %*% M %*% momentum
    energy[i,2] <- proposed_energy

    acceptance_ratio <- exp(current_energy - proposed_energy)
    metropolis_acceptance[i,1] <- stats::runif(1)
    proposals[i,] <- proposed_state

    if (metropolis_acceptance[i, 1] < acceptance_ratio) {
      current_state <- proposed_state
      samples[i,] <- current_state
      n_accepted <- n_accepted + 1
    }

    if (adaptive_mode) {
      if (acceptance_ratio > adaptive_target_acceptance) {
        step_size = step_size * 1.2
      } else {
        step_size = step_size * 1/1.2
      }
    }
  }

  return(list(samples=samples,
              proposals = proposals,
              energy = energy,
              metropolis_acceptance=metropolis_acceptance))
}
