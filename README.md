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
library(hamiltonianmcmc)
library(torch)
library(MASS)

d = 3
rho <- runif(3)
sigma <- c(3, 2, 1)
df <- mvrnorm(100, rep(0, d), t(rho%*%diag(sigma))%*%t(rho))

tgt_dens <- hamiltonianmcmc::TargetDensity$new(data = torch_tensor(df),
                                               init_mu = rnorm(d),
                                               diag(d))

num_samples <- 1500
initial_state <- rnorm(d)
step_size <- 0.0005
num_steps <- 30

hmc_res <- hamiltonian_mcmc(initial_state = initial_state, num_samples = num_samples,
                            step_size = step_size, num_steps = num_steps,
                            tgt_density = tgt_dens, M=rep(0.01, d),
                            adaptive_mode = T, adaptive_target_acceptance = 0.8)

plot(hmc_res$samples[100:1500,])
```
