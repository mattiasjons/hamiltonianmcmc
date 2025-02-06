library(R6)

#' @export
HMCAdapter <- R6Class("HMCAdapter",
                          public = list(
                            # Fields for count, mean, and M2 (for covariance calculations)
                            count = NULL,
                            metric_method = NULL,
                            metric_adapter = NULL,
                            metric_settings = NULL,
                            dual_adapter = NULL,
                            tmp_samples = NULL,
                            step_size = NULL,
                            reg_method = NULL,
                            reg_adapter = NULL,

                            # Initialization method
                            initialize = function(metric_method='welford', reg_method = 'none',
                                                  step_size, metric, metric_settings) {

                              # Initialize the Dual Averaging Adapter
                              self$step_size = step_size
                              self$dual_adapter = DualAveragingAdaptation$new(0.65, step_size)

                              self$metric_adapter <- DummyAdapter$new(metric)
                              self$count <- 0
                              self$metric_method = metric_method

                              self$metric_settings = metric_settings
                              self$reg_method = reg_method
                              self$reg_adapter = RegularizationAdapter$new(0.99, 20, 3000)
                            },

                            adapt_step = function(acceptance_prob) {
                              self$dual_adapter$adapt_step(acceptance_prob)
                            },

                            # Method to update the statistics with a new multidimensional value
                            add_sample = function(new_value, sjd=NULL) {

                              if(!self$metric_method =='welford') {

                                if (self$count <= self$metric_settings$k) {
                                  new_mt <- matrix(new_value, nrow = 1)
                                  self$tmp_samples <- rbind(self$tmp_samples, new_mt)
                                }

                                if (self$count == 2) {

                                  self$metric_adapter = WelfordAdapter$new(self$tmp_samples)

                                } else if (self$count > 2 & self$count < self$metric_settings$k) {

                                  self$metric_adapter$add_sample(new_value)

                                } else if (self$count==self$metric_settings$k) {

                                  self$metric_adapter <- switch(self$metric_method,
                                         "incpca" = IncPCA_Adapter$new(self$tmp_samples,
                                                                        self$metric_settings$k,
                                                                        3,
                                                                       self$reg_method),
                                         "ccipca" = CCIPCA_Adapter$new(self$tmp_samples,
                                                                       self$metric_settings$k,
                                                                       3,
                                                                       self$reg_method)
                                         )

                                  self$dual_adapter <- DualAveragingAdaptation$new(0.65,
                                                                                   self$step_size) #Reset dual averaging. Needed/useful?

                                } else if (self$count>self$metric_settings$k) {

                                  self$metric_adapter$add_sample(new_value)

                                  #Note: adapting squared jumping distance, DA
                                  self$reg_adapter$adapt_step_sjd(sjd)

                                  #Note: adapting using log likelihood
                                  #self$reg_adapter$adapt_step(self$get_eigvals(), self$get_eigvecs())
                                }
                              } else {

                                if (self$count <= 2) {
                                  new_mt <- matrix(new_value, nrow = 1)
                                  self$tmp_samples <- rbind(self$tmp_samples, new_mt)
                                }

                                if (self$count == 2) {

                                  self$metric_adapter = WelfordAdapter$new(self$tmp_samples)

                                } else if (self$count > 2) {

                                  self$metric_adapter$add_sample(new_value)
                                }
                              }

                              self$count <- self$count + 1
                            },

                            sample_covariance = function() {
                              if(self$metric_method %in% c('ccipca', 'incpca') && self$count>self$metric_settings$k) {
                                return(self$metric_adapter$sample_covariance(tau = self$reg_adapter$get_tau()))
                              } else {
                                return(self$metric_adapter$sample_covariance())
                              }
                            },

                            metric = function() {
                              if(self$metric_method %in% c('ccipca', 'incpca') && self$count>self$metric_settings$k) {
                                return(self$metric_adapter$metric(tau = self$reg_adapter$get_tau()))
                              } else {
                                return(self$metric_adapter$metric())
                              }
                            },

                            get_eigvals = function() {
                              if(self$metric_method %in% c('ccipca', 'incpca') && self$count>self$metric_settings$k) {
                                return(self$metric_adapter$get_eigvals())
                              } else {
                                cat('returned NA eigenvalues')
                                cat('\r\n')
                                return(rep(NA, self$metric_settings$k))
                              }
                            },

                            get_reg_eigvals = function() {
                              return(self$metric_adapter$get_reg_eigvals(tau = self$reg_adapter$get_tau()))
                            },

                            get_eigvecs = function() {
                              if(self$metric_method %in% c('ccipca', 'incpca') && self$count>self$metric_settings$k) {
                                return(self$metric_adapter$get_eigvecs())
                              } else {
                                cat('returned NA eigenvector')
                                cat('\r\n')
                                return(NA)
                              }
                            },

                            get_epsilon = function() {
                              return(self$dual_adapter$get_epsilon())
                            },

                            get_tau = function() {
                              return(self$reg_adapter$get_tau())
                            }
                          )
)

