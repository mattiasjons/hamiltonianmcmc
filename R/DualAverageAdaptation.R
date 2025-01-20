library(R6)

DualAveragingAdaptation <- R6Class("DualAveragingAdaptation",
                                   public = list(
                                     gamma = NULL,
                                     t0 = NULL,
                                     kappa = NULL,
                                     epsilon = NULL,
                                     mu_epsilon = NULL,
                                     log_epsilon_bar = NULL,
                                     h_bar_epsilon = NULL,
                                     adapt_epsilon_counter = NULL,
                                     target_accept_prob = NULL

                                     initialize = function(target_accept_prob, init_epsilon) {
                                       # Dual averaging constants (defaults similar to Stan)
                                       self$gamma <- 0.05
                                       self$t0 <- 10.0
                                       self$kappa <- 0.75

                                       # Initialize epsilon
                                       self$epsilon <- init_epsilon

                                       # Initialize mu values (centered logit_tau and logit_tau_2)
                                       self$mu_epsilon <- log(init_epsilon)

                                       # Initialize smoothed values
                                       self$log_epsilon_bar <- 0

                                       # Initialize h_bar values
                                       self$h_bar_epsilon <- 0.0

                                       # Initialize counters
                                       self$adapt_epsilon_counter <- 0

                                       # Target acceptance probability
                                       self$target_accept_prob <- target_accept_prob
                                     },

                                     adapt_step = function(accept_prob) {
                                       if (accept_prob > 1) {
                                         accept_prob <- 1
                                       }

                                       # Increment counters
                                       self$adapt_epsilon_counter <- self$adapt_epsilon_counter + 1

                                       # Compute eta for epsilon, logit_tau, and logit_tau_2
                                       eta_epsilon <- 1.0 / (self$adapt_epsilon_counter + self$t0)

                                       # Update h_bar for epsilon
                                       self$h_bar_epsilon <- (1 - eta_epsilon) * self$h_bar_epsilon + eta_epsilon * (self$target_accept_prob - accept_prob)

                                       # Update log epsilon
                                       log_epsilon <- self$mu_epsilon - (sqrt(self$adapt_epsilon_counter) / self$gamma) * self$h_bar_epsilon
                                       x_eta_epsilon <- self$adapt_epsilon_counter^(-self$kappa)
                                       self$log_epsilon_bar <- x_eta_epsilon * log_epsilon + (1 - x_eta_epsilon) * self$log_epsilon_bar
                                       self$epsilon <- exp(log_epsilon)
                                     },

                                     final_epsilon = function() {
                                       exp(self$log_epsilon_bar)
                                     },

                                     get_epsilon = function() {
                                       self$epsilon
                                     }
                                   ),
                                   private = list(
                                   )
)
