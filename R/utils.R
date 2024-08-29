#' Plot Leapfrog integration steps
#'
#' @param hmc_res A list containing the entries 'samples', 'leapfrog_states', normally output from the hamiltonian_mcmc function when returnleapfrogdetails=TRUE.
#' @param s Which specific sample to plot. If NULL, a random sample will be plotted.
#' @param d Which two variables will be plotted, given as a vector of two integers. If NULL, the two dimensions will be chosen at random.
#' @export
plot_leapfrog_steps <- function(hmc_res, s=NULL, d=NULL) {
  require(ggplot2)
  n <- nrow(hmc_res$samples)
  k <- ncol(hmc_res$samples)

  df <- MASS::mvrnorm(10000, colMeans(hmc_res$samples[(n/2):n,]), cov(hmc_res$samples[(n/2):n,])) #Posterior Predictive for Gaussian

  if (is.null(s)) {
    s <- sample((n/2):n, 1, F)
  }

  leapfrog_states <- do.call(rbind, hmc_res$leapfrog_states[[s]])

  if (is.null(d)) {
    smp_var <- sample(k, 2, F)
  } else {
    smp_var <- d
  }

  smp_leap <- leapfrog_states[,smp_var]

  ggplot(as.data.frame(df), aes(x=.data[[colnames(as.data.frame(df))[smp_var[1]]]], y=.data[[colnames(as.data.frame(df))[smp_var[2]]]]) ) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette= "Spectral", direction=1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      legend.position='none'
    ) +
    geom_path(mapping=aes(x=.data[[colnames(as.data.frame(leapfrog_states))[smp_var[1]]]], y=.data[[colnames(as.data.frame(leapfrog_states))[smp_var[2]]]]), data=as.data.frame(leapfrog_states), inherit.aes = F) +
    geom_point(mapping=aes(x=.data[[colnames(as.data.frame(leapfrog_states))[smp_var[1]]]], y=.data[[colnames(as.data.frame(leapfrog_states))[smp_var[2]]]]), data=as.data.frame(leapfrog_states)[1,], inherit.aes = F)
}

