library(R6)

DualAveragingAdaptation <- R6Class("DualAveragingAdaptation",
                                   public = list(
                                     initial_window = NULL,
                                     base_window = NULL,
                                     terminal_window = NULL,
                                     counter = NULL,
                                     gamma = NULL,
                                     t0 = NULL,
                                     kappa = NULL,
                                     epsilon = NULL,
                                     mu = NULL,
                                     log_epsilon_bar = NULL,
                                     h_bar = NULL,
                                     adapt_epsilon_counter = NULL,
                                     mass_matrix = NULL,
                                     inv_mass_matrix = NULL,
                                     target_accept_prob = NULL,
                                     warmup_steps = NULL,
                                     next_window = NULL,
                                     adapting = NULL,
                                     samples = NULL,
                                     num_samples = NULL,
                                     initialize = function(num_warmup_steps, target_accept_prob, init_epsilon, init_inv_mass_matrix) {
                                       # Windowing constants (defaults similar to Stan)
                                       self$initial_window = 75
                                       self$base_window <- 25
                                       self$terminal_window <- 50

                                       # Windowing counter
                                       self$counter <- 0

                                       # Dual averaging constants (defaults similar to Stan)
                                       self$gamma <- 0.05
                                       self$t0 <- 10.0
                                       self$kappa <- 0.75

                                       # Variables for dual averaging
                                       self$epsilon <- init_epsilon
                                       self$mu <- 0
                                       self$log_epsilon_bar <- 0
                                       self$h_bar <- 0.0
                                       self$adapt_epsilon_counter <- 0

                                       # Target acceptance probability
                                       self$target_accept_prob <- target_accept_prob

                                       # Warmup steps
                                       minimum_warmup_steps <- self$initial_window + self$terminal_window + self$base_window
                                       if (num_warmup_steps < minimum_warmup_steps) {
                                         stop(paste("Number of warmup steps less than the minimum value", minimum_warmup_steps))
                                       }

                                       self$warmup_steps <- num_warmup_steps
                                       self$next_window <- self$initial_window + self$base_window
                                       self$adapting <- TRUE
                                     },

                                     adapt_epsilon = function(accept_prob) {
                                       if (accept_prob > 1) {
                                         accept_prob <- 1
                                       }

                                       self$adapt_epsilon_counter <- self$adapt_epsilon_counter + 1

                                       eta <- 1.0 / (self$adapt_epsilon_counter + self$t0)

                                       self$h_bar <- (1 - eta) * self$h_bar + eta * (self$target_accept_prob - accept_prob)

                                       log_epsilon <- self$mu - (sqrt(self$adapt_epsilon_counter) / self$gamma) * self$h_bar

                                       x_eta <- self$adapt_epsilon_counter^(-self$kappa)
                                       self$log_epsilon_bar <- x_eta * log_epsilon + (1 - x_eta) * self$log_epsilon_bar
                                       self$epsilon <- exp(log_epsilon)
                                     },

                                     final_epsilon = function() {
                                       exp(self$log_epsilon_bar)
                                     },

                                     get_epsilon = function() {
                                       self$epsilon
                                     },

                                     step = function(x, accept_prob) {
                                       if (!self$adapting) {
                                         return(FALSE)
                                       }

                                       self$counter <- self$counter + 1

                                       if (self$counter >= self$warmup_steps) {
                                         self$epsilon <- self$final_epsilon()
                                         self$adapting <- FALSE
                                         return(FALSE)
                                       }

                                       self$adapt_epsilon(accept_prob)

                                       return(FALSE)
                                     }
                                   ),
                                   private = list(
                                   )
)
