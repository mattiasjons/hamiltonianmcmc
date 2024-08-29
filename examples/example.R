#install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(hamiltonianmcmc)
library(bayesplot)
library(cmdstanr)
color_scheme_set("brightblue")

#check_cmdstan_toolchain()
#install_cmdstan(cores = 2)

file <- file.path("./examples/normal.stan")
mod <- cmdstan_model(file)

data_list <- list(N = 10, y = rnorm(10, 0, 1))

hmc_res <- hamiltonian_mcmc(c(0, 1), 2000, 0.2, 10, stan_file = file,
                            stan_data = data_list, metric = diag(2),
                            metric_method = 'ccipca', adaptive_stepsize = T)

mcmc_intervals(hmc_res$samples[1000:2000,])
mcmc_trace(hmc_res$samples[500:2000,])
plot(log(hmc_res$step_sizes))

plot_leapfrog_steps(hmc_res)





load("./examples/data/mvnormal.rdata")

file <- file.path("./examples/mvnormal.stan")

k = 200
data_list <- list(N = nrow(X),
                  p = k,
                  y = y,
                  X = X[,order(sapply(1:ncol(X), function(i) cor(y, X[,i])), decreasing = T)[1:k]])

hmc_res <- hamiltonian_mcmc(rep(0, k), 500, 0.02, 20, stan_file = file, stan_data = data_list, metric = diag(k),
                            metric_method = 'ccipca', adaptive_stepsize = T)

mcmc_intervals(hmc_res$samples[200:500,])
mcmc_trace(hmc_res$samples[100:500,31:60])
plot(log(hmc_res$step_sizes))

plot(y, data_list$X %*% colMeans(hmc_res$samples[200:500,]))

plot_leapfrog_steps(hmc_res, s=5)
