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
                        l = NULL,
                        init_done = FALSE,
                        pca_values = NULL,
                        pca_vectors = NULL,
                        xbar = NULL,
                        variance_explained = NULL,
                        count = 0,

                        initialize = function(samples, variance_explained, l) {
                          require(onlinePCA)

                          stopifnot("Number of observations must be > 1" = nrow(samples)>1)

                          self$variance_explained <- variance_explained
                          self$l = l
                          pca <- eigen(cov(samples))
                          self$pca_values <- as.vector(pca$values)
                          self$pca_vectors <- pca$vectors
                          self$xbar <- colMeans(samples)
                          self$count <- nrow(samples)
                        },

                        # Method to update PCA with a new sample
                        add_sample = function(new_sample) {
                          new_sample <- as.vector(new_sample)
                          self$xbar <- self$update_mean(self$xbar, new_sample, self$count)

                          pca <- ccipca(self$pca_values, self$pca_vectors, new_sample, self$count,
                                        q = length(new_sample), center = self$xbar,
                                        l = min(self$l * 0.995^(self$count-1), self$count))

                          if (length(dim(pca$values))>1) {
                            values <- pca$values[,1]
                          } else {
                            values <- pca$values
                          }
                          self$pca_values <- values
                          self$pca_vectors <- pca$vectors

                          #self$pca_vectors <- pca$vectors #Needed?
                          self$count <- self$count + 1
                        },

                        sample_covariance = function() {
                          tmp <- which(cumsum(self$pca_values) / sum(self$pca_values) > self$variance_explained)[1]
                          lambda_shrunk <- ifelse(self$pca_values >= self$pca_values[tmp], self$pca_values, self$pca_values[tmp])

                          metric_inv <- self$pca_vectors %*% diag(lambda_shrunk) %*% t(self$pca_vectors)
                          return(metric_inv)
                        },

                        metric = function() {
                          tmp <- which(cumsum(self$pca_values) / sum(self$pca_values) > self$variance_explained)[1]
                          lambda_shrunk <- ifelse(self$pca_values >= self$pca_values[tmp], self$pca_values, self$pca_values[tmp])

                          metric <- self$pca_vectors %*% diag(1 / lambda_shrunk) %*% t(self$pca_vectors)
                          return(metric)
                        },

                        # Helper function to update the mean with a new sample
                        update_mean = function(xbar, new_sample, n) {
                          return(((n * xbar) + new_sample) / (n + 1))
                        }
                      ),
                      private = list(
                      )
)

