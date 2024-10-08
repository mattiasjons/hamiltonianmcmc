library(R6)

DualAveragingAdaptation <- R6Class("DualAveragingAdaptation",
                                   public = list(
                                     gamma = NULL,
                                     gamma_tau = NULL,
                                     t0 = NULL,
                                     t0_tau = NULL,
                                     kappa = NULL,
                                     epsilon = NULL,
                                     tau = NULL,              # tau in [0, 1]
                                     #tau_2 = NULL,            # tau_2 in [0, 1]
                                     logit_tau = NULL,        # logit-transformed tau
                                     #logit_tau_2 = NULL,      # logit-transformed tau_2
                                     mu_epsilon = NULL,
                                     mu_logit_tau = NULL,     # mu for logit tau
                                     #mu_logit_tau_2 = NULL,   # mu for logit tau_2
                                     log_epsilon_bar = NULL,
                                     logit_tau_bar = NULL,    # Smoothed logit_tau
                                     #logit_tau_2_bar = NULL,  # Smoothed logit_tau_2
                                     h_bar_epsilon = NULL,
                                     h_bar_logit_tau = NULL,  # h_bar for logit tau
                                     #h_bar_logit_tau_2 = NULL,  # h_bar for logit tau_2
                                     adapt_epsilon_counter = NULL,
                                     adapt_logit_tau_counter = NULL,   # Counter for logit tau adaptation
                                     #adapt_logit_tau_2_counter = NULL, # Counter for logit tau_2 adaptation
                                     target_accept_prob = NULL,

                                     initialize = function(target_accept_prob, init_epsilon, init_tau) {
                                       # Dual averaging constants (defaults similar to Stan)
                                       self$gamma <- 0.05
                                       self$t0 <- 10.0
                                       self$kappa <- 0.75

                                       # Initialize epsilon
                                       self$epsilon <- init_epsilon

                                       # Initialize tau and tau_2 (both in [0, 1])
                                       self$tau <- init_tau
                                       #self$tau_2 <- init_tau_2

                                       self$gamma_tau <- 0.15
                                       self$t0_tau <- 30.0

                                       # Initialize logit transformations of tau and tau_2
                                       self$logit_tau <- log(init_tau / (1 - init_tau))
                                       #self$logit_tau_2 <- log(init_tau_2 / (1 - init_tau_2))

                                       # Initialize mu values (centered logit_tau and logit_tau_2)
                                       self$mu_epsilon <- log(init_epsilon)
                                       self$mu_logit_tau <- self$logit_tau
                                       #self$mu_logit_tau_2 <- self$logit_tau_2

                                       # Initialize smoothed values
                                       self$log_epsilon_bar <- 0
                                       self$logit_tau_bar <- self$logit_tau
                                       #self$logit_tau_2_bar <- self$logit_tau_2

                                       # Initialize h_bar values
                                       self$h_bar_epsilon <- 0.0
                                       self$h_bar_logit_tau <- 0.0
                                       #self$h_bar_logit_tau_2 <- 0.0

                                       # Initialize counters
                                       self$adapt_epsilon_counter <- 0
                                       self$adapt_logit_tau_counter <- 0
                                       #self$adapt_logit_tau_2_counter <- 0

                                       # Target acceptance probability
                                       self$target_accept_prob <- target_accept_prob
                                     },

                                     adapt_step = function(accept_prob) {
                                       if (accept_prob > 1) {
                                         accept_prob <- 1
                                       }

                                       # Increment counters
                                       self$adapt_epsilon_counter <- self$adapt_epsilon_counter + 1
                                       self$adapt_logit_tau_counter <- self$adapt_logit_tau_counter + 1
                                       #self$adapt_logit_tau_2_counter <- self$adapt_logit_tau_2_counter + 1

                                       # Compute eta for epsilon, logit_tau, and logit_tau_2
                                       eta_epsilon <- 1.0 / (self$adapt_epsilon_counter + self$t0)
                                       eta_logit_tau <- 1.0 / (self$adapt_logit_tau_counter + self$t0_tau)
                                       #eta_logit_tau_2 <- 1.0 / (self$adapt_logit_tau_2_counter + self$t0_tau)

                                       # Update h_bar for epsilon
                                       self$h_bar_epsilon <- (1 - eta_epsilon) * self$h_bar_epsilon + eta_epsilon * (self$target_accept_prob - accept_prob)

                                       # Update log epsilon
                                       log_epsilon <- self$mu_epsilon - (sqrt(self$adapt_epsilon_counter) / self$gamma) * self$h_bar_epsilon
                                       x_eta_epsilon <- self$adapt_epsilon_counter^(-self$kappa)
                                       self$log_epsilon_bar <- x_eta_epsilon * log_epsilon + (1 - x_eta_epsilon) * self$log_epsilon_bar
                                       self$epsilon <- exp(log_epsilon)

                                       # Update h_bar for logit_tau. Add bias factor to
                                       self$h_bar_logit_tau <- (1 - eta_logit_tau) * self$h_bar_logit_tau + eta_logit_tau * (self$target_accept_prob - accept_prob)

                                       # Update logit_tau
                                       logit_tau <- self$mu_logit_tau - (sqrt(self$adapt_logit_tau_counter) / self$gamma_tau) * self$h_bar_logit_tau
                                       x_eta_logit_tau <- self$adapt_logit_tau_counter^(-self$kappa)
                                       self$logit_tau_bar <- x_eta_logit_tau * logit_tau + (1 - x_eta_logit_tau) * self$logit_tau_bar
                                       self$logit_tau <- logit_tau

                                       # Convert logit_tau back to [0, 1] using the sigmoid function
                                       self$tau <- 1 / (1 + exp(-self$logit_tau))

                                       # Update h_bar for logit_tau_2
                                       #self$h_bar_logit_tau_2 <- (1 - eta_logit_tau_2) * self$h_bar_logit_tau_2 + eta_logit_tau_2 * (self$target_accept_prob * 1.2 - accept_prob)

                                       # Update logit_tau_2
                                       #logit_tau_2 <- self$mu_logit_tau_2 - (sqrt(self$adapt_logit_tau_2_counter) / self$gamma_tau) * self$h_bar_logit_tau_2
                                       #x_eta_logit_tau_2 <- self$adapt_logit_tau_2_counter^(-self$kappa)
                                       #self$logit_tau_2_bar <- x_eta_logit_tau_2 * logit_tau_2 + (1 - x_eta_logit_tau_2) * self$logit_tau_2_bar
                                       #self$logit_tau_2 <- logit_tau_2

                                       # Convert logit_tau_2 back to [0, 1] using the sigmoid function
                                       #self$tau_2 <- 1 / (1 + exp(-self$logit_tau_2))
                                     },

                                     final_epsilon = function() {
                                       exp(self$log_epsilon_bar)
                                     },

                                     final_tau = function() {
                                       # Return final tau in [0, 1] by applying sigmoid to the smoothed logit_tau_bar
                                       1 / (1 + exp(-self$logit_tau_bar))
                                     },

                                     #final_tau_2 = function() {
                                      # # Return final tau_2 in [0, 1] by applying sigmoid to the smoothed logit_tau_2_bar
                                      # 1 / (1 + exp(-self$logit_tau_2_bar))
                                     #},

                                     get_epsilon = function() {
                                       self$epsilon
                                     },

                                     get_tau = function() {
                                       self$tau[1]  # tau is already transformed to be in [0, 1]
                                     }#,

                                     #get_tau_2 = function() {
                                    #   self$tau_2  # tau_2 is already transformed to be in [0, 1]
                                     #}
                                   ),
                                   private = list(
                                   )
)
