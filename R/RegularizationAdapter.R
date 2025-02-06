library(R6)

#' @export
RegularizationAdapter <- R6Class("RegularizationAdapter",
                                   public = list(

                                     tau = NA,
                                     min_f_tau = 1e30,
                                     init_temp = NULL,
                                     temp = NULL,
                                     n_samples = NULL,
                                     n_samples_seen = 0,
                                     lastk_f_tau_prop = NULL,
                                     steps_since_accept = NULL,
                                     adapt_sjd_counter = NULL,
                                     gamma = NULL,
                                     t0 = NULL,
                                     kappa = NULL,
                                     mu_tau = NULL,
                                     logit_tau_bar = NULL,
                                     h_bar_tau = NULL,
                                     esjd_short = NULL,
                                     esjd_long = NULL,

                                     initialize = function(tau = 0.99, init_temp = 25, n_samples) { #20 is good when constant
                                       self$min_f_tau <- 1e30
                                       self$tau <- tau
                                       self$init_temp <- init_temp
                                       self$temp <- init_temp
                                       self$n_samples <- n_samples
                                       self$lastk_f_tau_prop = rep(1e30, 50)
                                       self$steps_since_accept = 1
                                       self$adapt_sjd_counter = 0

                                       self$gamma <- 2
                                       self$t0 <- 100.0
                                       self$kappa <- 0.9

                                       # Initialize smoothed values
                                       self$logit_tau_bar <- 3
                                       self$mu_tau <- 3
                                       # Initialize h_bar values
                                       self$h_bar_tau <- 0.0

                                       self$esjd_short = rep(7, 3)
                                       self$esjd_long = rep(7, 50)
                                     },

                                     adapt_step = function(eigenvals, eigenvecs) {
                                       #browser()
                                       logit_tau <- log(self$tau/(1-self$tau))
                                       logit_tau_prop <- logit_tau + rnorm(1, mean = 0, sd = 0.3)
                                       tau_prop <- 1/(1+exp(-logit_tau_prop))
                                       f_tau_prop <- self$f_tau(eigenvals, eigenvecs, tau_prop)
                                       self$lastk_f_tau_prop <- c(f_tau_prop, self$lastk_f_tau_prop[1:49])

                                       f_tau_prop_star <- (0.8^self$steps_since_accept)*f_tau_prop + (1 - 0.8^self$steps_since_accept)*quantile(self$lastk_f_tau_prop, 0.2)

                                       cat(paste0('f_tau_prop: ', f_tau_prop, ', Min f_tau: ', self$min_f_tau, ', f_tau_star: ', f_tau_prop_star))
                                       cat('\r\n')

                                       #if (exp(-(f_tau_prop-self$min_f_tau)/self$temp)>runif(1)) {
                                       if (exp(-(f_tau_prop_star-self$min_f_tau)/self$temp)>runif(1)) {
                                         self$min_f_tau <- f_tau_prop_star
                                         self$tau <- tau_prop
                                         self$steps_since_accept = 0
                                         cat('Accepted new tau')
                                         cat('\r\n')
                                       }

                                       self$steps_since_accept = self$steps_since_accept + 1
                                       self$n_samples_seen <- self$n_samples_seen + 1
                                       #self$temp <- self$init_temp - self$n_samples_seen * ((self$init_temp - 0.01) / self$n_samples)
                                       self$temp <- self$init_temp * (0.001 / self$init_temp)^((self$n_samples_seen) / self$n_samples)
                                     },
                                     adapt_step_sjd = function(sjd) {
                                       if (is.na(sjd)) {
                                         sjd <- 0
                                       }

                                       # Increment counters
                                       self$adapt_sjd_counter <- self$adapt_sjd_counter + 1

                                       # Compute eta for epsilon, logit_tau, and logit_tau_2
                                       eta_tau <- 1.0 / (self$adapt_sjd_counter + self$t0)

                                       self$esjd_short <- c(sjd, self$esjd_short[1:2])
                                       cat('esjd_short: ')
                                       cat(mean(self$esjd_short))
                                       cat('\r\n')

                                       self$esjd_long <- c(sjd, self$esjd_long[1:49])
                                       cat('esjd_long: ')
                                       cat(quantile(self$esjd_long, 0.8))
                                       cat('\r\n')

                                       # Update h_bar for epsilon
                                       #self$h_bar_tau <- (1 - eta_tau) * self$h_bar_tau + eta_tau * (mean(self$esjd_short) - quantile(self$esjd_long, 0.8))
                                       self$h_bar_tau <- (1 - eta_tau) * self$h_bar_tau + eta_tau * (mean(self$esjd_short) - 7.5)

                                       # Update log epsilon
                                       logit_tau <- self$mu_tau - (sqrt(self$adapt_sjd_counter) / self$gamma) * self$h_bar_tau
                                       cat('logit_tau: ')
                                       cat(logit_tau)
                                       cat('\r\n')

                                       x_eta_tau <- self$adapt_sjd_counter^(-self$kappa)

                                       self$logit_tau_bar <- x_eta_tau * logit_tau + (1 - x_eta_tau) * self$logit_tau_bar
                                       self$tau <- 1 / (1 + exp(-logit_tau))

                                       cat('self_tau: ')
                                       cat(self$tau)
                                       cat('\r\n')
                                     },

                                     get_tau = function() {
                                       return(self$tau)
                                     },

                                     f_tau = function(eigenvals, eigenvecs, tau) {
                                       return(sum(diag(eigenvecs %*% diag(eigenvals) %*% t(eigenvecs) %*%
                                                  eigenvecs %*% diag(1/private$regularize_eigvals(eigenvals, tau)) %*% t(eigenvecs)))
                                       - log(prod(1/private$regularize_eigvals(eigenvals, tau))))
                                     }
                                   ),
                                   private = list(
                                     regularize_eigvals = function(eig_vals, tau) {
                                       if (tau==1) {
                                         return(eig_vals)
                                       } else {
                                         tmp <- which((cumsum(eig_vals) / sum(eig_vals)) > tau)[1]
                                         z_shrunk <- ifelse(eig_vals >= eig_vals[tmp], eig_vals, eig_vals[tmp])
                                         z_shrunk
                                       }
                                     }
                                   )
)
