library(R6)

MetricAdapter <- R6Class("MetricAdapter",
                                     public = list(),
                                     private = list()
)

#' @export
WelfordAdapter <- R6Class("WelfordAdapter",
                          inherit = MetricAdapter,
                          public = list(
                            # Fields for count, mean, and M2 (for covariance calculations)
                            count = NULL,
                            mean = NULL,
                            M2 = NULL,

                            # Initialization method
                            initialize = function(samples) {
                              stopifnot("Number of observations must be > 1" = nrow(samples)>1)

                              self$mean <- colMeans(samples)
                              self$count <- nrow(samples)

                              deviations <- sweep(samples, 2, self$mean)

                              M2 <- t(deviations) %*% deviations
                              self$M2 <- M2
                            },

                            # Method to update the statistics with a new multidimensional value
                            add_sample = function(new_value) {
                              stopifnot('Adapter must have been initialized with a dataset with >1 observation' = self$count>1,
                                        'Observation must be numeric' = is.numeric(new_value))
                              self$count <- self$count + 1

                              new_value <- as.vector(new_value)

                              delta <- new_value - self$mean
                              self$mean <- self$mean + delta / self$count
                              delta2 <- new_value - self$mean
                              self$M2 <- self$M2 + outer(delta, delta2)
                            },

                            # Method to return sample variance
                            sample_covariance = function() {
                              stopifnot('Adapter must have been initialized with a dataset with >1 observation' = self$count>1)
                              sample_covariance <- self$M2 / (self$count - 1)
                              sample_covariance <- (self$count / (self$count + 5.0)) * sample_covariance +
                                1e-3 * (5.0 / (self$count + 5.0)) * diag(nrow(sample_covariance))
                              return(sample_covariance)
                            },

                            metric = function() {
                              stopifnot('Adapter must have been initialized with a dataset with >1 observation' = self$count>1)
                              covar <- self$sample_covariance()
                              metric_inv <- (self$count / (self$count + 5.0)) * covar +
                                1e-3 * (5.0 / (self$count + 5.0)) * diag(nrow(covar))
                              metric <- solve(metric_inv)
                              return(metric)
                            }
                          )
                          )

#' @export
CCIPCA_Adapter <- R6Class("CCIPCA_Adapter",
                          inherit = MetricAdapter,
                      public = list(
                        k = NULL,
                        l = NULL,
                        init_done = FALSE,
                        pca_values = NULL,
                        pca_vectors = NULL,
                        xbar = NULL,
                        count = 0,
                        reg_method = NULL,

                        initialize = function(samples, k, l, reg_method='none') {
                          require(onlinePCA)

                          stopifnot("Number of observations must be > 1" = nrow(samples)>1)

                          self$k = k
                          self$l = l
                          pca <- eigen(cov(samples))
                          self$pca_values <- as.vector(pca$values)[1:k]
                          self$pca_vectors <- pca$vectors[,1:k]
                          self$xbar <- colMeans(samples)
                          self$count <- nrow(samples)
                          self$reg_method <- reg_method
                        },

                        # Method to update PCA with a new sample
                        add_sample = function(new_sample) {
                          new_sample <- as.vector(new_sample)
                          self$xbar <- self$update_mean(self$xbar, new_sample, self$count)

                          pca <- ccipca(self$pca_values, self$pca_vectors, new_sample, self$count,
                                        q = self$k, center = self$xbar,
                                        #l = min(self$l * 0.995^(self$count-1), self$count))
                                        #l = min(4, self$count))
                                        l = 0)

                          if (length(dim(pca$values))>1) {
                            values <- pca$values[,1]
                          } else {
                            values <- pca$values
                          }

                          self$pca_values <- values
                          self$pca_vectors <- pca$vectors
                          self$count <- self$count + 1
                        },

                        sample_covariance = function(beta = 0.5, tau_min = 0.01, tau_max = 0.99) {
                          lambda_shrunk <- NA
                          switch (self$reg_method,
                            'none' = {
                              lambda_shrunk <- self$pca_values
                              },
                            'linear' = {
                              lambda_shrunk <- (1 - beta) * self$pca_values + beta * mean(self$pca_values)
                              },
                            'minmax' = {

                              start_eigval <- ceiling(tau_min * length(self$pca_values))

                              lambda_shrunk <- self$pca_values
                              lambda_shrunk[1:start_eigval] <- lambda_shrunk[start_eigval]

                              tmp <- which((cumsum(lambda_shrunk) / sum(lambda_shrunk)) > tau_max)[1]
                              lambda_shrunk <- ifelse(lambda_shrunk >= lambda_shrunk[tmp], lambda_shrunk, lambda_shrunk[tmp])

                              },
                            'truncate' = {
                              lambda_shrunk <- ifelse(self$pca_values > 0, self$pca_values, self$pca_values[self$pca_values>0])
                              }
                          )

                          metric_inv <- self$pca_vectors %*% diag(lambda_shrunk) %*% t(self$pca_vectors)
                          return(metric_inv)
                        },

                        metric = function(beta = 0.5, tau_min = 0.01, tau_max = 0.99) {
                          lambda_shrunk <- NA
                          switch (self$reg_method,
                                  'none' = {
                                    lambda_shrunk <- self$pca_values
                                  },
                                  'linear' = {
                                    lambda_shrunk <- (1 - beta) * self$pca_values + beta * mean(self$pca_values)
                                  },
                                  'minmax' = {

                                    start_eigval <- ceiling(tau_min * length(self$pca_values))

                                    lambda_shrunk <- self$pca_values
                                    lambda_shrunk[1:start_eigval] <- lambda_shrunk[start_eigval]

                                    tmp <- which((cumsum(lambda_shrunk) / sum(lambda_shrunk)) > tau_max)[1]
                                    lambda_shrunk <- ifelse(lambda_shrunk >= lambda_shrunk[tmp], lambda_shrunk, lambda_shrunk[tmp])

                                  },
                                  'truncate' = {
                                    lambda_shrunk <- ifelse(self$pca_values > 0, self$pca_values, self$pca_values[self$pca_values>0])
                                  }
                          )

                          metric <- self$pca_vectors %*% diag(1 / lambda_shrunk) %*% t(self$pca_vectors)
                          return(metric)
                        },

                        # Helper function to update the mean with a new sample
                        update_mean = function(xbar, new_sample, n) {
                          return(((n * xbar) + new_sample) / (n + 1))
                        },

                        get_eigvals = function() {
                          return(self$pca_values)
                        },

                        get_eigvecs = function() {
                          return(self$pca_vectors)
                        }
                      ),
                      private = list(
                      )
)

#' @export
IncPCA_Adapter <- R6Class("IncPCA_Adapter",
                          inherit = MetricAdapter,
                          public = list(
                            k = NULL, #Number of components
                            l = NULL, #Amnesic factor
                            init_done = FALSE,
                            pca_values = NULL,
                            pca_vectors = NULL,
                            xbar = NULL,
                            count = 0,
                            reg_method = NULL,

                            initialize = function(samples, k, l, reg_method='none') {
                              require(onlinePCA)

                              stopifnot("Number of observations must be > 1" = nrow(samples)>1)

                              self$k = k
                              self$l = l
                              pca <- eigen(cov(samples))
                              self$pca_values <- as.vector(pca$values)[1:k]
                              self$pca_vectors <- pca$vectors[,1:k]
                              self$xbar <- colMeans(samples)
                              self$count <- nrow(samples)
                              self$reg_method <- reg_method
                            },

                            # Method to update PCA with a new sample
                            add_sample = function(new_sample) {
                              #incRpca(lambda, U, x, n, f = 1/n, q = length(lambda), center, tol = 1e-7)
                              new_sample <- as.vector(new_sample)
                              self$xbar <- self$update_mean(self$xbar, new_sample, self$count)

                              pca <- incRpca(lambda=self$pca_values, U=self$pca_vectors,
                                             x=new_sample, n=self$count,
                                             q = self$k, center = self$xbar,
                                             f = 0)

                              if (length(dim(pca$values))>1) {
                                values <- pca$values[,1]
                              } else {
                                values <- pca$values
                              }

                              self$pca_values <- values
                              self$pca_vectors <- pca$vectors
                              self$count <- self$count + 1
                            },

                            sample_covariance = function(beta = 0.5, tau_min = 0.01, tau_max = 0.99) {
                              lambda_shrunk <- NA
                              switch (self$reg_method,
                                      'none' = {
                                        lambda_shrunk <- self$pca_values
                                      },
                                      'linear' = {
                                        lambda_shrunk <- (1 - beta) * self$pca_values + beta * mean(self$pca_values)
                                      },
                                      'minmax' = {
                                        start_eigval <- ceiling(tau_min * length(self$pca_values))

                                        lambda_shrunk <- self$pca_values
                                        lambda_shrunk[1:start_eigval] <- lambda_shrunk[start_eigval]

                                        tmp <- which((cumsum(lambda_shrunk) / sum(lambda_shrunk)) > tau_max)[1]
                                        lambda_shrunk <- ifelse(lambda_shrunk >= lambda_shrunk[tmp], lambda_shrunk, lambda_shrunk[tmp])

                                      },
                                      'truncate' = {
                                        lambda_shrunk <- ifelse(self$pca_values > 0, self$pca_values, self$pca_values[self$pca_values>0])
                                      }
                              )

                              metric_inv <- self$pca_vectors %*% diag(lambda_shrunk) %*% t(self$pca_vectors)
                              return(metric_inv)
                            },

                            metric = function(beta = 0.5, tau_min = 0.01, tau_max = 0.99) {
                              lambda_shrunk <- NA
                              switch (self$reg_method,
                                      'none' = {
                                        lambda_shrunk <- self$pca_values
                                      },
                                      'linear' = {
                                        lambda_shrunk <- (1 - beta) * self$pca_values + beta * mean(self$pca_values)
                                      },
                                      'minmax' = {

                                        start_eigval <- ceiling(tau_min * length(self$pca_values))

                                        lambda_shrunk <- self$pca_values
                                        lambda_shrunk[1:start_eigval] <- lambda_shrunk[start_eigval]

                                        tmp <- which((cumsum(lambda_shrunk) / sum(lambda_shrunk)) > tau_max)[1]
                                        lambda_shrunk <- ifelse(lambda_shrunk >= lambda_shrunk[tmp], lambda_shrunk, lambda_shrunk[tmp])

                                      },
                                      'truncate' = {
                                        lambda_shrunk <- ifelse(self$pca_values > 0, self$pca_values, self$pca_values[self$pca_values>0])
                                      }
                              )

                              metric <- self$pca_vectors %*% diag(1 / lambda_shrunk) %*% t(self$pca_vectors)
                              return(metric)
                            },

                            # Helper function to update the mean with a new sample
                            update_mean = function(xbar, new_sample, n) {
                              return(((n * xbar) + new_sample) / (n + 1))
                            },

                            get_eigvals = function() {
                              return(self$pca_values)
                            },

                            get_eigvecs = function() {
                              return(self$pca_vectors)
                            }
                          ),
                          private = list(
                          )
)

#' @export
DummyAdapter <- R6Class("DummyAdapter",
                          inherit = MetricAdapter,
                          public = list(
                            metric_var = NULL,

                            # Initialization method
                            initialize = function(metric) {
                              self$metric_var = metric
                            },

                            # Method to update the statistics with a new multidimensional value
                            add_sample = function(new_value) {
                            },

                            # Method to return sample variance
                            sample_covariance = function() {
                              return(solve(self$metric_var))
                            },

                            metric = function() {
                              return(self$metric_var)
                            }
                          )
)
