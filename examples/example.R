#install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
library(hamiltonianmcmc)
library(bayesplot)
library(cmdstanr)
color_scheme_set("brightblue")

#check_cmdstan_toolchain()
#install_cmdstan(cores = 2)


eg_fn <- function() rexp(1)
df <- generate_eig_df(301, 0.9, 25, eg_fn)
X <- df[[1]]
y <- df[[2]]

file <- file.path("./examples/mvnormal.stan")


k = 300
data_list <- list(N = nrow(X),
                  p = k,
                  y = y,
                  X = X[,order(sapply(1:ncol(X), function(i) cor(y, X[,i])), decreasing = T)[1:k]])

library(glmnet) #Let's "cheat" to find good starting values.

#RegularizationAdapter$debug('adapt_step')
#RegularizationAdapter$undebug('adapt_step')

hmc_res <- hamiltonian_mcmc(c(as.numeric(coef.glmnet(glmnet(data_list$X, y, intercept = F, lambda = 0.01)))[-1], 1),
                            num_samples = 2000, 0.1, 20, stan_file = file, stan_data = data_list, metric = diag(k+1),
                            metric_method = 'ccipca', metric_adapter_settings = list(k=50, reg_method='minmax'))

hmc_res_inc <- hamiltonian_mcmc(c(as.numeric(coef.glmnet(glmnet(data_list$X, y, intercept = F, lambda = 0.01)))[-1], 1),
                            num_samples = 2000, 0.1, 20, stan_file = file, stan_data = data_list, metric = diag(k+1),
                            metric_method = 'incpca', metric_adapter_settings = list(k=50, reg_method='minmax'))

hmc_res_cov <- hamiltonian_mcmc(c(as.numeric(coef.glmnet(glmnet(data_list$X, y, intercept = F, lambda = 0.01)))[-1], 1),
                            2000, 0.1, 20, stan_file = file, stan_data = data_list, metric = diag(k+1),
                            metric_method = 'welford')

mcmc_intervals(hmc_res$samples[500:1000,])
mcmc_intervals(hmc_res_inc$samples[500:1000,])
mcmc_intervals(hmc_res_cov$samples[500:1000,])

mcmc_trace(hmc_res$samples[1000:2000,1:30])
acf(hmc_res$samples[1000:2000,115])

mcmc_trace(hmc_res_inc$samples[1000:2000,1:30])
acf(hmc_res_inc$samples[1000:2000,115])

mcmc_trace(hmc_res_cov$samples[300:500,1:30])
acf(hmc_res_cov$samples[1000:2000,115])
#acf(hmc_res_inc$samples[seq(1, 500, by=2),80])

hist(rowSums(diff(hmc_res_inc$samples)^2))
hist(rowSums(diff(hmc_res_inc$samples[seq(2, 500, by=2),])^2))
mcmc_trace(hmc_res_cov$samples[1000:2000,1:30])

hmc_lst <- list(Spectral=hmc_res, Spectral_Incremental=hmc_res_inc, Welford=hmc_res_cov)
plot_step_size(hmc_lst)
plot_posterior_draws(hmc_lst, 1000)

plot_eigenvalues_minmax_shrink(hmc_res, T)
plot_eigenvalues_minmax_shrink(hmc_res_inc, T)

plot_tau(hmc_res)
plot_tau(hmc_res_inc)

plot_leapfrog_steps(hmc_res)
plot_leapfrog_steps(hmc_res_inc)
plot_leapfrog_steps(hmc_res_cov)

mean(hmc_res$ess)
mean(hmc_res_cov$ess)
mean(hmc_res_inc$ess)

plot(hmc_res$tau[-1,1], sqrt(get_esjd(hmc_res)[,2]))

plot_ess_var(hmc_lst)
plot_ess(hmc_lst)
plot_esjd(hmc_lst)

plot_trace(hmc_lst)

plot(get_esjd(hmc_res), ylim=c(0, 1000))
esjd_ts <- lowess(get_esjd(hmc_res)[,2], f = 0.05)
lines(esjd_ts$y, col='red', lwd=3)

#plot(get_condition_number(hmc_res, T))
#plot(get_condition_number(hmc_res_inc, T))
plot_esjd_potential(hmc_res)
plot_esjd_trace(hmc_lst)

esjd_ts_cov <- lowess(get_esjd(hmc_res_cov)[,2], f = 0.05)
lines(esjd_ts_cov$y, col='red', lwd=3)

eff_smp <- zoo::rollapply(hmc_res$samples, 50, coda::effectiveSize, by=10)
eff_smp_inc <- zoo::rollapply(hmc_res_inc$samples, 50, coda::effectiveSize, by=10)
eff_smp_cov <- zoo::rollapply(hmc_res_cov$samples, 50, coda::effectiveSize, by=10)

df1 <- data.frame(e=rowMeans(eff_smp), t=approx(hmc_res$ess_s[,1], n = 196)$y, Method='Spectral')
df2 <- data.frame(e=rowMeans(eff_smp_cov), t=approx(hmc_res_cov$ess_s[,1], n = 196)$y, Method='Welford')
df3 <- data.frame(e=rowMeans(eff_smp_inc), t=approx(hmc_res_inc$ess_s[,1], n = 196)$y, Method='Spectral Incremental')
df1$t <- df1$t-min(df1$t)
df2$t <- df2$t-min(df2$t)
df3$t <- df3$t-min(df3$t)

df_plot <- rbind(df1, df2, df3)
ggplot(data=df_plot, mapping = aes(x=t, y=e, col=Method, group=Method)) +
  geom_line() + geom_smooth() + labs(x='Time (s)', y='Effective Sample Size, rolling estimate (n=30)')




reg_eig <- get_reg_eig(hmc_res)
eig_vec <- hmc_res$eig_vectors[[131]]

plot(eig_vec %*% diag(1/reg_eig[reg_eig$k==131,]$ev) %*% t(eig_vec), hmc_res$mass_matrices[[131]])

#View(hmc_res$eig_vectors[[401]] - hmc_res$eig_vectors[[300]])
plot(eig_vec %*% diag(1/hmc_res$eig_values[131,]) %*% t(eig_vec), hmc_res$mass_matrices[[132]])

tau_rng <- runif(10000, 0.98, 1)
hist(tau_rng)
tau_rng <- tau_rng[sample(1:10000, 2000, F)]
tau_rng <- tau_rng[tau_rng<1]

eig_vec <- hmc_res$eig_vectors[[500]]

esjd_tau <- sapply(tau_rng, function(tau_max_1) {
  smp_pt <- sample(451:500, 1)
  cat(tau_max_1)
  cat('\r\n')
  M <- eig_vec %*% diag(1/regularize_eigvals(hmc_res$eig_values[500,], 0.01, tau_max_1)) %*% t(eig_vec)
  M_inv <- eig_vec %*% diag(regularize_eigvals(hmc_res$eig_values[500,], 0.01, tau_max_1)) %*% t(eig_vec)

  p <- MASS::mvrnorm(1, rep(0, 500), M_inv)

  potential_energy <- -as.numeric(hmc_res$stan_obj$log_prob(c(hmc_res$samples[smp_pt,])))
  kinetic_energy <- 0.5 * t(p) %*% M_inv %*% p
  current_energy <- potential_energy + kinetic_energy

  lf <- leapfrog_integration(c(hmc_res$samples[smp_pt,]),
                             momentum = p,
                             step_size = hmc_res$step_sizes[500],
                             stan_obj = hmc_res$stan_obj,
                             num_steps = 20, returndetails = F, M)

  proposed_potential_energy <- -as.numeric(hmc_res$stan_obj$log_prob(lf$proposed_state))
  proposed_kinetic_energy <- 0.5 * t(lf$proposed_momentum) %*% M_inv %*% lf$proposed_momentum
  proposed_energy <- proposed_potential_energy + proposed_kinetic_energy

  acceptance_ratio <- min(exp(current_energy - proposed_energy), 1)

  acceptance_ratio * sum((lf$proposed_state-c(hmc_res$samples[smp_pts[1],]))^2)
})

plot(tau_rng, esjd_tau)

esjd_slices <- sapply(seq(0.981, 0.999, 0.0002), function(x) {
  esjd_slice <- esjd_tau[tau_rng>(x-0.0001) & tau_rng<(x+0.0001)]
  length(esjd_slice[esjd_slice<1e-32])/length(esjd_slice)
})

esjd_var_slices <- sapply(seq(0.981, 0.999, 0.0002), function(x) {
  esjd_slice <- esjd_tau[tau_rng>(x-0.0001) & tau_rng<(x+0.0001)]
  var(esjd_slice)
})

qplot(y=esjd_var_slices, x=1:91) + geom_point() + geom_smooth(span=0.75)

size_slices <- sapply(seq(0.981, 0.999, 0.0002), function(x) {
  esjd_slice <- esjd_tau[tau_rng>(x-0.0001) & tau_rng<(x+0.0001)]
  length(esjd_slice)
})

plot(x=seq(0.981, 0.999, 0.0002), y=esjd_slices, pch=16)

