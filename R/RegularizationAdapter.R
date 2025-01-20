library(R6)

RegularizationAdapter <- R6Class("RegularizationAdapter",
                                   public = list(

                                     tau_min = NA,
                                     tau_max = NA,
                                     min_ll = 1e30,

                                     initialize = function(tau_min=0.1, tau_max=0.9) {
                                       self$min_ll <- Inf
                                       self$tau_min <- tau_min
                                       self$tau_max <- tau_max
                                     },

                                     adapt_step = function(eigenvals, eigenvecs) {
                                       c(sum(diag(eigenvecs %*% diag(self$regularize_eigvals(eigenvals, tau[1], tau[2])) %*% t(hmc_res$eig_vectors[[i]]) %*% hmc_res$eig_vectors[[i]] %*% diag(1/regularize_eigvals(hmc_res$eig_values[i,], tau[1], tau[2])) %*% t(hmc_res$eig_vectors[[i]]))),
                                         -log(prod(1/regularize_eigvals(hmc_res$eig_values[i,], tau[1], tau[2]))))
                                     }
                                   ),
                                   private = list(
                                     regularize_eigvals = function(eig_vals, tau_min, tau_max) {

                                         start_eigval <- ceiling(tau_min * length(eig_vals))
                                         z_shrunk <- eig_vals
                                         z_shrunk[1:start_eigval] <- z_shrunk[start_eigval]

                                         tmp <- which((cumsum(z_shrunk) / sum(z_shrunk)) > tau_max)[1]
                                         z_shrunk <- ifelse(z_shrunk >= z_shrunk[tmp], z_shrunk, z_shrunk[tmp])
                                         z_shrunk
                                     }
                                   )
)
