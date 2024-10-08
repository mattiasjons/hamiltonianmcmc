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




eg_fn <- function() rexp(1)
df <- generate_eig_df(501, 0.8, 25, eg_fn)
X <- df[[1]]
y <- df[[2]]

file <- file.path("./examples/mvnormal.stan")


k = 500
data_list <- list(N = nrow(X),
                  p = k,
                  y = y,
                  X = X[,order(sapply(1:ncol(X), function(i) cor(y, X[,i])), decreasing = T)[1:k]])

library(glmnet) #Let's "cheat" to find good starting values.

hmc_res <- hamiltonian_mcmc(as.numeric(coef.glmnet(glmnet(data_list$X, y, intercept = F, lambda = 0.01)))[-1],
                            num_samples = 500, 0.1, 20, stan_file = file, stan_data = data_list, metric = diag(k),
                            metric_method = 'ccipca', metric_adapter_settings = list(k=84))

hmc_res_cov <- hamiltonian_mcmc(as.numeric(coef.glmnet(glmnet(data_list$X, y, intercept = F, lambda = 0.01)))[-1],
                            500, 0.1, 20, stan_file = file, stan_data = data_list, metric = diag(k),
                            metric_method = 'welford')

mcmc_intervals(hmc_res$samples[300:500,])
mcmc_intervals(hmc_res_cov$samples[300:500,])

mcmc_trace(hmc_res$samples[300:500,1:30])
mcmc_trace(hmc_res_cov$samples[300:500,1:30])

hmc_lst <- list(Spectral=hmc_res, Welford=hmc_res_cov)
plot_step_size(hmc_lst)

plot_eigenvalues_linear_shrink(hmc_res, T)

plot_tau(hmc_res)

plot_leapfrog_steps(hmc_res_cov)
plot_leapfrog_steps(hmc_res)

mean(hmc_res$ess)
mean(hmc_res_cov$ess)

plot_ess(hmc_lst)
plot_esjd(hmc_lst)
