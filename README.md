# hamiltonianmcmc

R package to run Hamiltonian MCMC sampling

# How to Install

Easiest way to install is:

```r
devtools::install_github("mattiasjons/hamiltonianmcmc", upgrade_dependencies = FALSE)
```

# How to use

A simple example for how to use the code follows below:

```r
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
df <- generate_eig_df(501, 0.9, 25, eg_fn)
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
                            metric_method = 'ccipca', metric_adapter_settings = list(k=50, reg_method='minmax'))

hmc_res_inc <- hamiltonian_mcmc(as.numeric(coef.glmnet(glmnet(data_list$X, y, intercept = F, lambda = 0.01)))[-1],
                            num_samples = 500, 0.1, 20, stan_file = file, stan_data = data_list, metric = diag(k),
                            metric_method = 'incpca', metric_adapter_settings = list(k=50, reg_method='minmax'))

hmc_res_cov <- hamiltonian_mcmc(as.numeric(coef.glmnet(glmnet(data_list$X, y, intercept = F, lambda = 0.01)))[-1],
                            500, 0.1, 20, stan_file = file, stan_data = data_list, metric = diag(k),
                            metric_method = 'welford')

mcmc_intervals(hmc_res$samples[300:500,])
mcmc_intervals(hmc_res_inc$samples[300:500,])
mcmc_intervals(hmc_res_cov$samples[300:500,])

mcmc_trace(hmc_res$samples[300:500,1:30])
mcmc_trace(hmc_res_inc$samples[300:500,1:30])
mcmc_trace(hmc_res_cov$samples[300:500,1:30])

plot(density(hmc_res$samples[400:500, 4]))
lines(density(hmc_res_cov$samples[400:500, 4]))
lines(density(hmc_res_inc$samples[400:500, 4]))

hmc_lst <- list(Spectral=hmc_res, Spectral_Incremental=hmc_res_inc, Welford=hmc_res_cov)
plot_step_size(hmc_lst)

plot_eigenvalues_minmax_shrink(hmc_res, T)
plot_eigenvalues_minmax_shrink(hmc_res_inc, T)


plot_tau(hmc_res)
plot(hmc_res$tau[,2])
plot(-log(1/hmc_res_inc$tau[,2]-1), hmc_res_inc$step_sizes)


plot_leapfrog_steps(hmc_res_cov)
plot_leapfrog_steps(hmc_res)

mean(hmc_res$ess)
mean(hmc_res_cov$ess)
mean(hmc_res_inc$ess)

plot_ess_var(hmc_lst)
plot_ess(hmc_lst)
plot_esjd(hmc_lst)


esjd <- data.frame(i = 1:(nrow(hmc_res$samples)-1), esjd=rowSums(diff(hmc_res$samples)^2))
plot(hmc_res$energy[-1,3], esjd$esjd)

plot(get_condition_number(hmc_res, T))

```
