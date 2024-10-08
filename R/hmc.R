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
#' @param stan_obj A number.
#' @param num_steps A number.
#' @export
#' @returns A list of proposed state, proposed momentum and final gradient
leapfrog_integration <- function(current_state, momentum,
                                 step_size, stan_obj,
                                 num_steps, returndetails=F, M) {

  states = vector("list", num_steps + 1)
  momenta = vector("list", num_steps + 2)

  proposed_momentum <- momentum + 0.5 * c(step_size) * as.numeric(stan_obj$grad_log_prob(current_state))
  proposed_state <- current_state + c(step_size) * M %*% proposed_momentum

  states[[1]] <- c(current_state)
  momenta[[1]] <- c(proposed_momentum)

  for (j in 1:num_steps) {
    if(any(is.na(proposed_state)) | any(is.na(proposed_momentum))) {
      return(list(diverge = T))
    }
    proposed_momentum <- proposed_momentum + c(step_size) * as.numeric(stan_obj$grad_log_prob(proposed_state))
    proposed_state <- proposed_state + c(step_size) * M %*% proposed_momentum

    states[[j + 1]] <- c(proposed_state)
    momenta[[j + 1]] <- c(proposed_momentum)
  }

  final_grad <- as.numeric(stan_obj$grad_log_prob(proposed_state))
  proposed_momentum <- proposed_momentum + 0.5 * c(step_size) * final_grad
  momenta[[num_steps + 2]] <- c(proposed_momentum)

  return(list(diverge = F, proposed_state = proposed_state, proposed_momentum = proposed_momentum,
              gradient = final_grad, states = states, momenta = momenta))
}

#' Sample using hamiltonian MCMC
#'
#' @param initial_state A vector with length \eqn{d}, representing the initial state.
#' @param num_samples Number of samples to generate.
#' @param step_size How large is the step size in the leapfrof integration.
#' @param num_steps Number of leap frog steps. If Adaptive Mode is enabled, this represents the initial value.
#' @param stan_obj An object of class TargetDensity representing the target density to sample from.
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
                             stan_file,
                             stan_data,
                             metric_method='ccipca',
                             metric = diag(length(initial_state)),
                             adaptive_stepsize_target_acceptance = 0.65,
                             returnleapfrogdetails = F,
                             metric_adapter_settings = list(NA)) {
  #TODO: Asserts

  library(cmdstanr)
  require(posterior)
  require(coda)

  file <- file.path(stan_file)
  mod <- cmdstan_model(file, compile = F)
  mod$compile(compile_model_methods=TRUE, force_recompile=T)

  fit <- mod$sample(
    iter_warmup = 1,
    iter_sampling = 1,
    data = stan_data,
    #adapt_delta = 0.99,
    seed = 123,
    chains = 4,
    parallel_chains = 4, max_treedepth = 14,
    refresh = 0,
    save_warmup = 1,
    metric='dense_e',
    save_metric = T,
    show_messages = F,
    show_exceptions = F,
    diagnostics = ''
  )

  cat('Model compiled and prepared\r\n')

  if(is.null(metric)) {
    metric <- rep(1, length(initial_state))
  }

  if (length(metric)==length(initial_state) & length(metric)>1) {
    metric = diag(metric)
  }

  metric_inv <- MASS::ginv(metric)

  eg_val <- NA
  if (metric_method =='ccipca') {
    library(onlinePCA)
    eg_val <- matrix(nrow = num_samples, ncol=metric_adapter_settings$k)
  }

  current_state <- initial_state
  samples <- matrix(nrow = num_samples, ncol = length(initial_state))
  rejected <- matrix(nrow = 100 * num_samples, ncol = length(initial_state)) #TODO: Make this more efficient...
  energy <- matrix(nrow = num_samples, ncol=2)
  leapfrogmomenta <- list()
  leapfrogstates <- list()
  metropolis_acceptance = matrix(nrow = num_samples, ncol=1)
  step_sizes = matrix(nrow = num_samples, ncol=1)
  mass_matrices <- list()
  init_done <- F
  n_accepted <- 0
  last_updated_mass_step <- 0
  ess_s <- matrix(nrow = num_samples/5, ncol=2)
  tau_hist <- matrix(nrow = num_samples, ncol=2)

  #adapter <- DualAveragingAdaptation$new(0.65, step_size, 0.2)
  adapter <- HMCAdapter$new(metric_method=metric_method,
                        step_size=step_size, metric=metric,
                        metric_settings = metric_adapter_settings)

  i = 1
  n_rejected = 1

  while (i <= num_samples) {
    if (i%%50==0) {
      cat(paste0('Generated ', i, ' samples.'))
      cat('\r\n')
    }

    accepted = F
    n_inner_iter = 0
    while(!accepted) {
      n_inner_iter = n_inner_iter + 1
      momentum <- MASS::mvrnorm(1, mu = rep(0, length(current_state)), Sigma = metric)

      # Leapfrog integration
      leapfrog_result <- leapfrog_integration(current_state, momentum, adapter$get_epsilon(), fit, num_steps, returnleapfrogdetails, metric_inv)
      if (leapfrog_result$diverge) {
        proposed_state <- NA
        proposed_momentum <- NA
      } else {
        proposed_state <- leapfrog_result$proposed_state
        proposed_momentum <- leapfrog_result$proposed_momentum
      }

      # Metropolis acceptance step
      if(any(is.na(proposed_state)) | any(is.na(proposed_momentum))) {
        current_energy <- -as.numeric(fit$log_prob(current_state)) + 0.5 * t(momentum) %*% metric_inv %*% momentum
        proposed_energy <- NA
        acceptance_ratio <- 0
      } else {
        current_energy <- -as.numeric(fit$log_prob(current_state)) + 0.5 * t(momentum) %*% metric_inv %*% momentum
        proposed_energy <- -as.numeric(fit$log_prob(proposed_state)) + 0.5 * t(proposed_momentum) %*% metric_inv %*% proposed_momentum

        acceptance_ratio <- exp(current_energy - proposed_energy)
        if (is.na(acceptance_ratio)) {
          acceptance_ratio <- 0
        }
        if (is.infinite(acceptance_ratio)) {
          acceptance_ratio <- 1
        }
      }

      adapter$adapt_step(acceptance_ratio)
      cat(adapter$get_epsilon())
      cat('\r\n')

      U <- stats::runif(1)

      if (!is.na(proposed_energy) && !is.na(acceptance_ratio) && U < acceptance_ratio) {

        adapter$add_sample(proposed_state)
        metric <- adapter$metric()
        metric_inv <- adapter$sample_covariance()

        current_state <- proposed_state
        samples[i,] <- current_state
        leapfrogmomenta[[i]] <- leapfrog_result$momenta
        leapfrogstates[[i]] <- leapfrog_result$states
        energy[i,1] <- current_energy
        energy[i,2] <- proposed_energy
        metropolis_acceptance[i,1] <- U
        mass_matrices[[i]] <- metric

        accepted = T

        if (i%%5==0) {
          ess_s[i/5,] <- c(Sys.time(), mean(apply(samples[(i-4):i,], 2, effectiveSize)))
        }

      } else {
        rejected[n_rejected,] <- proposed_state
        n_rejected <- n_rejected + 1
        if (n_rejected%%50==0) {
          cat(paste0('Rejected ', n_rejected, ' samples.'))
          cat('\r\n')
        }
      }
    }

    step_sizes[i] = adapter$get_epsilon()

    if (metric_method =='ccipca') {
      tau_hist[i, ] <- c(adapter$get_tau(), adapter$get_tau())
      eg_val[i,] <- adapter$get_eigvals()
    }

    i = i + 1

  }

  ess <- effectiveSize(samples)

  return(list(samples=as_draws_matrix(samples),
              energy = energy,
              metropolis_acceptance=metropolis_acceptance,
              mass_matrices = mass_matrices,
              step_sizes = step_sizes,
              leapfrog_momenta = leapfrogmomenta,
              leapfrog_states = leapfrogstates,
              rejected = rejected, #TODO: return only divergent samples instead. All rejected are probably not interesting.
              stan_obj = fit,
              ess = ess,
              ess_s = ess_s,
              tau = tau_hist,
              eig_values = eg_val
              ))
}


