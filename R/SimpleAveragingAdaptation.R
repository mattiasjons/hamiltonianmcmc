library(R6)

SimpleAveragingAdaptation <- R6Class("SimpleAveragingAdaptation",
                                   public = list(
                                     adapt_epsilon_counter = NULL,
                                     mass_matrix = NULL,
                                     inv_mass_matrix = NULL,
                                     target_accept_prob = NULL,
                                     warmup_steps = NULL,
                                     next_window = NULL,
                                     adapting = NULL,
                                     t0 = NULL,
                                     epsilon = NULL,
                                     h_bar = NULL,
                                     initialize = function(target_accept_prob, init_epsilon) {
                                       # Windowing counter
                                       self$adapt_epsilon_counter <- 0

                                       # Dual averaging constants (defaults similar to Stan)
                                       #self$gamma <- 0.05
                                       self$t0 <- 10.0
                                       self$h_bar <- 0
                                       #self$kappa <- 0.75

                                       # Variables for dual averaging
                                       self$epsilon <- init_epsilon

                                       # Target acceptance probability
                                       self$target_accept_prob <- target_accept_prob

                                       self$adapting <- TRUE
                                     },

                                     adapt_epsilon = function(accept_prob) {
                                       if (accept_prob > 1) {
                                         accept_prob <- 1
                                       }

                                       self$adapt_epsilon_counter <- self$adapt_epsilon_counter + 1

                                       eta <- 0.99

                                       self$h_bar <- (1 - eta) * self$h_bar + eta * (accept_prob - self$target_accept_prob)

                                       self$epsilon <- self$epsilon + self$epsilon * sign(self$h_bar) * 0.01
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
