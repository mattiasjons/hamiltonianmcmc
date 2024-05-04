library(R6)
library(torch)

#' @export
TargetDensity <- R6Class("TargetDensity",
                      public = list(
                        data = NULL,
                        mu = NULL,
                        sigma = NULL,
                        initialize = function(data, init_mu, init_sigma) {
                          self$data = data
                          self$mu = torch_tensor(init_mu, requires_grad = T)
                          self$sigma = torch_tensor(init_sigma)
                        },
                        value = function(mu) {
                          self$mu$set_data(mu)
                          mvn = torch::distr_multivariate_normal(
                            loc = self$mu,
                            covariance_matrix = self$sigma
                          )
                          log_likelihood = torch_sum(mvn$log_prob(self$data))
                          return(log_likelihood)
                        },
                        grad = function(mu) {
                          ll = self$value(mu)
                          return(autograd_grad(ll, self$mu))
                        }))
