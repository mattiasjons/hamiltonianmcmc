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

k = 500
data_list <- list(N = nrow(X),
                  p = k,
                  y = y,
                  X = X[,order(sapply(1:ncol(X), function(i) cor(y, X[,i])), decreasing = T)[1:k]])

library(glmnet) #Let's "cheat" to find good starting values.

hmc_res <- hamiltonian_mcmc(as.numeric(coef.glmnet(glmnet(data_list$X, y, intercept = F, lambda = 1)))[-1],
                            500, 0.1, 20, stan_file = file, stan_data = data_list, metric = diag(k),
                            metric_method = 'ccipca', adaptive_stepsize = T, 4)

hmc_res_cov <- hamiltonian_mcmc(as.numeric(coef.glmnet(glmnet(data_list$X, y, intercept = F, lambda = 1)))[-1],
                            500, 0.1, 20, stan_file = file, stan_data = data_list, metric = diag(k),
                            metric_method = 'cov', adaptive_stepsize = T)

mcmc_intervals(hmc_res$samples[300:500,])
mcmc_intervals(hmc_res_cov$samples[300:500,])

mcmc_trace(hmc_res$samples[300:500,1:30])
mcmc_trace(hmc_res_cov$samples[300:500,1:30])

step_size_eigen <- log(hmc_res$step_sizes)
step_size_cov <- log(hmc_res_cov$step_sizes)
plot(1:500, step_size_eigen, col=1, pch=16,
     xlab='Sample', ylab='log(Step size)')
points(1:500, step_size_cov, col=2, pch=16)
legend('bottomright', legend=c('Welford', 'Spectral'), col=c('red', 'black'), pch=16)

plot(y, data_list$X %*% colMeans(hmc_res$samples[200:500,]))

plot_leapfrog_steps(hmc_res_cov, s=30)

mean(hmc_res$ess)
mean(hmc_res_cov$ess)

ess_s <- hmc_res$ess_s
ess_s[,1] <- ess_s[,1] - min(ess_s[,1])

ess_s_cov <- hmc_res_cov$ess_s
ess_s_cov[,1] <- ess_s_cov[,1] - min(ess_s_cov[,1])

plot(ess_s[,1], cumsum(ess_s[,2])/ess_s[,1], type = 'l',
     lwd=2, xlim = c(0, max(ess_s_cov[,1])),
     ylim=c(0, 12), xlab='Seconds', ylab='ESS/s')

lines(ess_s_cov[,1],
      cumsum(ess_s_cov[,2])/ess_s_cov[,1],
      col='red', lwd=2)

