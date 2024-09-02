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


calculate_ess <- function(samples) {
  n <- length(samples)

  autocorrelation <- calculate_autocorrelation_fft(samples)

  sum_autocorr <- 0
  for (k in 1:(n - 1)) {
    if (autocorrelation[k + 1] < 0) {
      break
    }
    sum_autocorr <- sum_autocorr + autocorrelation[k + 1]
  }

  ess <- n / (1 + 2 * sum_autocorr)

  return(ess)
}

calculate_autocorrelation_fft <- function(x) {
  n <- length(x)

  # Subtract the mean
  x_centered <- x - mean(x)

  # FFT of the time series
  fft_x <- fft(x_centered)

  # Compute the power spectrum
  power_spectrum <- fft_x * Conj(fft_x)

  # Inverse FFT to get the autocovariance function
  autocovariance <- Re(fft(power_spectrum, inverse = TRUE) / n)

  # Normalize by the variance to get the autocorrelation function
  autocorrelation <- autocovariance / autocovariance[1]

  # Return the autocorrelation function
  return(autocorrelation)
}


adapt_epsilon <- function(accept_prob, adapt_epsilon_counter,
                          h_bar, target_accept_prob,
                          mu, gamma, kappa, t0,
                          log_epsilon_bar, epsilon) {

  # Ensure accept_prob does not exceed 1
  if (accept_prob > 1) {
    accept_prob <- 1.0
  }

  # Increment the adaptation counter
  adapt_epsilon_counter <- adapt_epsilon_counter + 1
  counter <- adapt_epsilon_counter

  # Calculate eta
  eta <- 1.0 / (counter + t0)

  # Update h_bar
  h_bar <- (1 - eta) * h_bar + eta * (target_accept_prob - accept_prob)

  # Calculate log_epsilon
  log_epsilon <- mu - (sqrt(counter) / gamma) * h_bar

  # Calculate x_eta
  x_eta <- counter^(-kappa)

  # Update log_epsilon_bar
  log_epsilon_bar <- x_eta * log_epsilon + (1 - x_eta) * log_epsilon_bar

  # Update epsilon
  epsilon <- exp(log_epsilon)

  return(list(
    epsilon = epsilon,
    h_bar = h_bar,
    log_epsilon_bar = log_epsilon_bar,
    adapt_epsilon_counter = adapt_epsilon_counter
  ))
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
                             adaptive_stepsize = TRUE,
                             adaptive_stepsize_target_acceptance = 0.65,
                             returnleapfrogdetails = F,
                             epsilon = 1e-2) {
  #TODO: Asserts

  library(cmdstanr)
  require(posterior)

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

  if (metric_method =='ccipca') {
    library(onlinePCA)
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

  #adapter <- SimpleAveragingAdaptation$new(target_accept_prob = adaptive_stepsize_target_acceptance,
  #                                       init_epsilon = step_size)
  adapter <- DualAveragingAdaptation$new(2000, 0.65, step_size, diag(length(initial_state)))

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

      if (adaptive_stepsize) {
        adapter$adapt_epsilon(acceptance_ratio)
      }

      U <- stats::runif(1)

      if (!is.na(proposed_energy) && !is.na(acceptance_ratio) && U < acceptance_ratio) {

        current_state <- proposed_state
        samples[i,] <- current_state
        leapfrogmomenta[[i]] <- leapfrog_result$momenta
        leapfrogstates[[i]] <- leapfrog_result$states
        energy[i,1] <- current_energy
        energy[i,2] <- proposed_energy
        metropolis_acceptance[i,1] <- U
        mass_matrices[[i]] <- metric

        if(metric_method=='ccipca'){
          if (!init_done && i > 1) {
            pca <- eigen(cov(samples[!is.na(samples[,1]),]))
            xbar <- colMeans(samples[!is.na(samples[,1]),])

            #pca <- list(values=pca$sdev^2,
            #            vectors=pca$loadings)

            init_done <- T
            last_updated_mass_step <- i

            #values <- length(ncol(samples)) * (pca$values + 1e-3) / sum(pca$values + 1e-3)
            values <- pca$values

            # Step 2: Shrink the eigenvalues
            #alpha <- 0.9^i  # Shrinkage parameter (0 <= alpha <= 1)
            #lambda_target <- rep(1, length(initial_state))  # Target value for shrinkage (e.g., mean of eigenvalues)

            # Apply the shrinkage
            #lambda_shrunk <- (i / (i + 5.0)) * values + 1e-3 * (5.0 / (i + 5.0))

            tmp <- which(cumsum(values)/sum(values) > 0.9999)[1]
            lambda_shrunk <- ifelse(values >= values[tmp], values, values[tmp])

            #metric_inv <- pca$vectors %*% diag(values) %*% t(pca$vectors)
            #metric <- pca$vectors %*% diag(1/values) %*% t(pca$vectors)
            metric_inv <- pca$vectors %*% diag(lambda_shrunk) %*% t(pca$vectors)
            metric <- pca$vectors %*% diag(1/lambda_shrunk) %*% t(pca$vectors)

          } else if (init_done) {
            smp <- samples[i,]

            xbar <- updateMean(xbar, smp, i-1)

            pca <- ccipca(pca$values, pca$vectors, smp, i-1, q = length(smp), center = xbar, l = min(i-1, 3))
            last_updated_mass_step <- i
            if (length(dim(pca$values))>1) {
              values <- pca$values[,1]
            } else {
              values <- pca$values
            }

            # Step 2: Shrink the eigenvalues
            #alpha <- 0.9^i  # Shrinkage parameter (0 <= alpha <= 1)
            #lambda_target <- rep(1, length(values))  # Target value for shrinkage (e.g., mean of eigenvalues)

            # Apply the shrinkage
            #lambda_shrunk <- (1 - alpha) * values + alpha * lambda_target
            #lambda_shrunk <- (i / (i + 5.0)) * values + 1e-3 * (5.0 / (i + 5.0))
            tmp <- which(cumsum(values)/sum(values) > 0.9999)[1]
            lambda_shrunk <- ifelse(values >= values[tmp], values, values[tmp])

            #metric_inv <- pca$vectors %*% diag(values) %*% t(pca$vectors)
            #metric <- pca$vectors %*% diag(1/values) %*% t(pca$vectors)
            metric_inv <- pca$vectors %*% diag(lambda_shrunk) %*% t(pca$vectors)
            metric <- pca$vectors %*% diag(1/lambda_shrunk) %*% t(pca$vectors)

            #values <- length(ncol(samples)) * (values + 1e-3) / sum(values + 1e-3)
            #metric_inv <- pca$vectors %*% diag(values) %*% t(pca$vectors)
            #metric <- pca$vectors %*% diag(1/(values)) %*% t(pca$vectors)
          }
        } else if (metric_method == 'cov' & i > 1) {
          covar <- cov(samples[!is.na(samples[,1]),])

          metric_inv <- (i / (i + 5.0)) * covar +
            1e-3 * (5.0 / (i + 5.0)) * diag(nrow(covar))
          metric <- solve(metric_inv)
        }

        accepted = T
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
    i = i + 1
  }

  ess <- apply(samples[!is.na(samples[,1]),], 2, calculate_ess)

  return(list(samples=as_draws_matrix(samples),
              energy = energy,
              metropolis_acceptance=metropolis_acceptance,
              mass_matrices = mass_matrices,
              step_sizes = step_sizes,
              leapfrog_momenta = leapfrogmomenta,
              leapfrog_states = leapfrogstates,
              rejected = rejected,
              ess = ess))
}


