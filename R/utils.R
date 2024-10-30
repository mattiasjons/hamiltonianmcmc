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
    stat_density_2d(aes(fill = after_stat(density)), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette= "Spectral", direction=1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
      legend.position='none'
    ) +
    geom_path(mapping=aes(x=.data[[colnames(as.data.frame(leapfrog_states))[smp_var[1]]]], y=.data[[colnames(as.data.frame(leapfrog_states))[smp_var[2]]]]), data=as.data.frame(leapfrog_states), inherit.aes = F) +
    geom_point(mapping=aes(x=.data[[colnames(as.data.frame(leapfrog_states))[smp_var[1]]]], y=.data[[colnames(as.data.frame(leapfrog_states))[smp_var[2]]]]), data=as.data.frame(leapfrog_states)[1,], inherit.aes = F)
}

#' Plot eigenvalue spectra for each step of the adaptation
#'
#' @param hmc_res A list containing the entry 'eig_values' normally output from the hamiltonian_mcmc function when metric_method is set to something that returns eigenvalues.
#' @param log Should the plot be drawn on log scale?
#' @param margin Should a marginal density of the lower eigenvalue treshold be included?
#' @export
plot_eigenvalues <- function(hmc_res, log=F, margin=F) {
  require(ggplot2)

  get_df <- function(k, log=F) {
    n <- ncol(hmc_res$eig_values)
    z_raw <- hmc_res$eig_values[k,]
    z_cs <- cumsum(z_raw)/sum(z_raw)
    z_0 <- which(z_cs>hmc_res$tau[k,1])[1]
    z_1 <- max(which(z_cs>hmc_res$tau[k,2])[1] - 1, 1)
    if (log) {
      z_df <- data.frame(k=k, i=1:n, ev=log(z_raw), trunc=1:n>z_0)
      z_cut <- data.frame(k=k, i=z_0:n, ev=log(z_raw[z_0]))
    } else {
      z_df <- data.frame(k=k, i=1:n, ev=z_raw, trunc=1:n>z_0)
      z_cut <- data.frame(k=k, i=z_0:n, ev=z_raw[z_0])
    }

    return(list(z_df, z_cut))
  }
  n <- nrow(hmc_res$samples)
  n_min <- which(!is.na(hmc_res$eig_values[,1]))[1]
  dfs <- lapply(n_min:n, get_df, log=log)
  z_df <- do.call(rbind, lapply(dfs, function(d) d[[1]]))
  z_cut <- do.call(rbind, lapply(dfs, function(d) d[[2]]))

  gg <- ggplot(z_df, aes(x=i, y=ev, col=trunc, group=k)) + geom_line(alpha=0.3) +
               geom_line(inherit.aes = F, data = z_cut, aes(x=i, y=ev, group=k), col='red', alpha=0.3) +
               scale_color_discrete(guide='none') +
               geom_point(inherit.aes = F, data=z_cut, aes(x=0, y=ev), alpha=0)
  if(margin) {
    require(ggExtra)
    ggMarginal(gg, margins = "y", col = "red", fill = "red", size = 8)
  } else {
    gg
  }
}

#' Plot eigenvalue spectra for each step of the adaptation. Use this alternative when linear shrinking has been applied
#'
#' @param hmc_res A list containing the entry 'eig_values' normally output from the hamiltonian_mcmc function when metric_method is set to something that returns eigenvalues.
#' @param log Should the plot be drawn on log scale?
#' @export
plot_eigenvalues_linear_shrink <- function(hmc_res, log=F) {
  require(ggplot2)

  get_df <- function(k, log=F) {
    n <- ncol(hmc_res$eig_values)
    z_raw <- hmc_res$eig_values[k,]
    z_shrunk <- (1-hmc_res$tau[k,]) * hmc_res$eig_values[k,] + hmc_res$tau[k,] * mean(hmc_res$eig_values[k,])

    if (log) {
      z_df <- data.frame(k=k, i=1:n, ev=log(z_raw))
      z_df_shrunk <- data.frame(k=k, i=1:n, ev=log(z_shrunk))
    } else {
      z_df <- data.frame(k=k, i=1:n, ev=z_raw)
      z_df_shrunk <- data.frame(k=k, i=1:n, ev=z_shrunk)
    }

    return(list(z_df, z_df_shrunk))
  }
  n <- nrow(hmc_res$samples)
  n_min <- which(!is.na(hmc_res$eig_values[,1]))[1]
  dfs <- lapply(n_min:n, get_df, log=log)
  z_df <- do.call(rbind, lapply(dfs, function(d) d[[1]]))
  z_df$Type <- 'Raw'
  z_shrunk <- do.call(rbind, lapply(dfs, function(d) d[[2]]))
  z_shrunk$Type <- 'Shrunk'

  z_df <- rbind(z_df, z_shrunk)

  ggplot(z_df, aes(x=i, y=ev, col=Type, group=interaction(k,Type))) + geom_line(alpha=0.3) +
    scale_color_discrete(guide='none')
}

#' Plot eigenvalue spectra for each step of the adaptation.
#' Use this alternative when minimum/maximum truncation/shrinking of eigenvalues has been applied.
#'
#' @param hmc_res A list containing the entry 'eig_values' normally output from the hamiltonian_mcmc function when metric_method is set to something that returns eigenvalues.
#' @param log Should the plot be drawn on log scale?
#' @export
plot_eigenvalues_minmax_shrink <- function(hmc_res, log=F) {
  require(ggplot2)

  get_df <- function(k, log=F) {
    n <- ncol(hmc_res$eig_values)
    z_raw <- hmc_res$eig_values[k,]
    start_eigval <- ceiling(hmc_res$tau[k,1] * n)

    z_shrunk <- z_raw
    z_shrunk[1:start_eigval] <- z_shrunk[start_eigval]

    tmp <- which((cumsum(z_shrunk) / sum(z_shrunk)) > hmc_res$tau[k,2])[1]
    z_shrunk <- ifelse(z_shrunk >= z_shrunk[tmp], z_shrunk, z_shrunk[tmp])

    if (log) {
      z_df <- data.frame(k=k, i=1:n, ev=log(z_raw))
      z_df_shrunk <- data.frame(k=k, i=1:n, ev=log(z_shrunk))
    } else {
      z_df <- data.frame(k=k, i=1:n, ev=z_raw)
      z_df_shrunk <- data.frame(k=k, i=1:n, ev=z_shrunk)
    }

    return(list(z_df, z_df_shrunk))
  }
  n <- nrow(hmc_res$samples)
  n_min <- which(!is.na(hmc_res$eig_values[,1]))[1]
  dfs <- lapply(n_min:n, get_df, log=log)
  z_df <- do.call(rbind, lapply(dfs, function(d) d[[1]]))
  z_df$Type <- 'Raw'
  z_shrunk <- do.call(rbind, lapply(dfs, function(d) d[[2]]))
  z_shrunk$Type <- 'Shrunk'

  z_df <- rbind(z_df, z_shrunk)

  ggplot(z_df, aes(x=i, y=ev, col=Type, group=interaction(k,Type))) + geom_line(alpha=0.3) +
    scale_color_discrete(guide='none')
}

#' Plot the adaptation of parameters tau and tau_2
#'
#' @param hmc_res A list containing the entry 'tau' normally output from the hamiltonian_mcmc function when metric_method is set to something that returns eigenvalues.
#' @export
plot_tau <- function(hmc_res){
  require(reshape2)
  require(ggplot2)

  tau_df <- melt(hmc_res$tau)
  tau_df$Variable <- ifelse(tau_df$Var2==1, 'tau[2]', 'tau')

  ggplot(tau_df, mapping = aes(x=Var1, y=value, group=Var2, col=Variable)) +
    scale_colour_discrete(name = "Variable", labels = expression(tau, tau[2])) +
    geom_point() +
    geom_smooth(span = 0.1)
}


#' Generate a sample dataset from a set of random eigenvalues
#'
#' @param k Total number of columns generated, i.e. number of columns of X + 1. Default=501
#' @param beta Expected proportion of eigenvalues that are 0. Default=0.7
#' @param nobs Number of observations in generated dataset
#' @param fn Function to generate a single eigenvalue
#' @export
generate_eig_df <- function(k=501, beta=0.7, nobs=25, fn){
  des <- sapply(1:k, function(z) ifelse(runif(1)>beta, fn(), 0))
  n <- length(des)
  s <- diag(des)

  q <- qr.Q(qr(matrix(runif(n*n), nrow=n)))
  semidef <- t(q) %*% s %*% q

  df <- MASS::mvrnorm(nobs, mu = rep(0, k), Sigma = semidef)
  X <- df[,1:(k-1)]
  y <- df[,k]
  return(list(X, y))
}

#' Plot Effective Sample Size / second
#'
#' @param hmc_res_list A list of hmc result objects
#' @export
plot_ess <- function(hmc_res_list) {
  require(ggplot2)

  hmc_lst_2 <- lapply(hmc_res_list, function(hmc_res) {
    ess_s <- hmc_res$ess_s
    ess_s[,1] <- ess_s[,1] - min(ess_s[,1])
    ess_s <- cbind(ess_s, cumsum(ess_s[,2])/ess_s[,1])
    ess_s
  })

  ess_df <- as.data.frame(do.call(rbind, hmc_lst_2))
  ess_df$method <- rep(names(hmc_lst_2), sapply(hmc_lst_2, nrow))

  ggplot(ess_df, aes(x=V1, y=V3, group=method, col=method)) +
    scale_x_continuous(name='Seconds') +
    scale_y_continuous(name='ESS/s') +
    scale_color_discrete(name='Method') +
    geom_line(size=1.2)
}


#' Plot Effective Sample Size per variable
#'
#' @param hmc_res_list A list of hmc result objects
#' @export
plot_ess_var <- function(hmc_res_list) {
  require(ggplot2)

  lst_ess <- lapply(hmc_lst, function(hmc_res) {
    data.frame(ess=hmc_res$ess)
  })

  df_ess <- as.data.frame(do.call(rbind, lst_ess))
  df_ess$Method <- rep(names(lst_ess), sapply(lst_ess, nrow))

  ggplot(df_ess, aes(x=ess, col=Method, fill=Method)) + geom_density(alpha=0.5)
}


#' Plot Step Size
#'
#' @param hmc_res_list A list of hmc result objects
#' @export
plot_step_size <- function(hmc_res_list) {
  require(ggplot2)

  hmc_lst_2 <- lapply(hmc_lst, function(hmc_res) {
    data.frame(s=1:length(hmc_res$step_sizes), epsilon=log(hmc_res$step_sizes))
  })

  step_size_df <- as.data.frame(do.call(rbind, hmc_lst_2))
  step_size_df$method <- rep(names(hmc_lst_2), sapply(hmc_lst_2, nrow))

  ggplot(step_size_df, aes(x=s, y=epsilon, group=method, col=method)) +
    scale_x_continuous(name='Step') +
    scale_y_continuous(name='Step size') +
    scale_color_discrete(name='Method') +
    geom_point()
}

#' Plot Expected Squared Jumping Distance
#'
#' @param hmc_res_list A list of hmc result objects
#' @export
plot_esjd <- function(hmc_lst) {
  esjd <- lapply(hmc_lst, function(hmc_res) {
    mh_cond <- ifelse(exp(rowSums(hmc_res$energy[,3:4])-rowSums(hmc_res$energy[,1:2]))>1,
                      1,
                      exp(rowSums(hmc_res$energy[,3:4])-rowSums(hmc_res$energy[,1:2])))
    data.frame(i = 1:(nrow(hmc_res$samples)-1), esjd=mh_cond[-1] * rowSums(diff(hmc_res$samples)^2))
  })
  esjd <- do.call(rbind, esjd)
  esjd$method <- rep(names(hmc_lst), sapply(hmc_lst, function(hmc_res) nrow(hmc_res$samples)-1))

  ggplot(esjd, aes(x=i, y=esjd, col=method)) + geom_point(alpha=0.2) + geom_smooth() +
    coord_cartesian(ylim = c(0, quantile(esjd$esjd, 0.99)))
}


#' Calculate the condition number, either based on the raw eigenvalues, or based on the regularized ones.
#' Note: Currently the regularized method only supports the 'minmax' regularization method.
#'
#' @param hmc_res A hmc result objects containing eig_values and tau entries
#' @param regularized Should the condition value be calculated based on regularized eigenvalues
#' @param reg_method If regularized eigenvalues are used, which method has been used for regularization
#' @export
get_condition_number <- function(hmc_res, regularized=F, reg_method='minmax') {

  #TODO: Assert that the hmc list contains eig_values, and if regularized that reg_method matches
  if(regularized && reg_method=='minmax') {
    get_df <- function(k) {
      n <- ncol(hmc_res$eig_values)
      z_raw <- hmc_res$eig_values[k,]
      start_eigval <- ceiling(hmc_res$tau[k,1] * n)

      z_shrunk <- z_raw
      z_shrunk[1:start_eigval] <- z_shrunk[start_eigval]

      tmp <- which((cumsum(z_shrunk) / sum(z_shrunk)) > hmc_res$tau[k,2])[1]
      z_shrunk <- ifelse(z_shrunk >= z_shrunk[tmp], z_shrunk, z_shrunk[tmp])

      z_df_shrunk <- data.frame(k=k, i=1:n, ev=z_shrunk)

      return(z_df_shrunk)
    }

    n <- nrow(hmc_res$samples)
    n_min <- which(!is.na(hmc_res$eig_values[,1]))[1]
    dfs <- lapply(n_min:n, get_df)
    z_df <- do.call(rbind, dfs)
    return(aggregate(z_df$ev, by=list(z_df$k), FUN=max)$x/
             aggregate(z_df$ev, by=list(z_df$k), FUN=min)$x)
  } else {
    apply(hmc_res$eig_values[!is.na(hmc_res$eig_values[,1]),], 1, max)/
      apply(hmc_res$eig_values[!is.na(hmc_res$eig_values[,1]),], 1, min)
  }
}



#' Get regularized eigenvalues
#' Note: Currently the regularized method only supports the 'minmax' regularization method.
#'
#' @param hmc_res A hmc result objects containing eig_values and tau entries.
#' @param reg_method If regularized eigenvalues are used, which method has been used for regularization
#' @export
get_reg_eig <- function(hmc_res, reg_method='minmax') {

  #TODO: Assert that the hmc list contains eig_values, and if regularized that reg_method matches
  if(reg_method=='minmax') {
    get_df <- function(k) {
      n <- ncol(hmc_res$eig_values)
      z_raw <- hmc_res$eig_values[k,]
      start_eigval <- ceiling(hmc_res$tau[k,1] * n)

      z_shrunk <- z_raw
      z_shrunk[1:start_eigval] <- z_shrunk[start_eigval]

      tmp <- which((cumsum(z_shrunk) / sum(z_shrunk)) > hmc_res$tau[k,2])[1]
      z_shrunk <- ifelse(z_shrunk >= z_shrunk[tmp], z_shrunk, z_shrunk[tmp])

      z_df_shrunk <- data.frame(k=k, i=1:n, ev=z_shrunk)

      return(z_df_shrunk)
    }

    n <- nrow(hmc_res$samples)
    n_min <- which(!is.na(hmc_res$eig_values[,1]))[1]
    dfs <- lapply(n_min:n, get_df)
    z_df <- do.call(rbind, dfs)
    return(z_df)
  } else {
    return(NA)
  }
}


#' Calculate Squared Jumping Distance
#'
#' @param hmc_res A hmc result objects containing energy and samples entries.
#' @export
get_esjd <- function(hmc_res) {
  mh_cond <- ifelse(exp(rowSums(hmc_res$energy[,3:4])-rowSums(hmc_res$energy[,1:2]))>1,
                    1,
                    exp(rowSums(hmc_res$energy[,3:4])-rowSums(hmc_res$energy[,1:2])))
  data.frame(i = 1:(nrow(hmc_res$samples)-1), esjd=mh_cond[-1] * rowSums(diff(hmc_res$samples)^2))
}


#' Plot Trace (sum of diagonal/eigenvalues of sample covariance matrix)
#'
#' @param hmc_lst A list of hmc result objects.
#' @export
plot_trace <- function(hmc_lst) {

  traces <- lapply(hmc_lst, function(hmc_res) {
    if (length(hmc_res$eig_values)==1) { #If length=1, eig_values contains a single NA entry.
      trace_increment <- unlist(lapply(hmc_res$mass_matrices, function(M) sum(diag(solve(M)))))
      data.frame(k = 1:length(trace_increment), trace=trace_increment)
    } else {
      reg_eig_df <- get_reg_eig(hmc_res, 'minmax')
      reg_eig_df <- aggregate(reg_eig_df$ev, by=list(reg_eig_df$k), sum)
      colnames(reg_eig_df) <- c('k', 'trace')
      reg_eig_df
    }
  })

  trace_df <- do.call(rbind, traces)
  rownames(trace_df) <- NULL

  trace_df$method <- rep(names(hmc_lst), unlist(lapply(traces, nrow)))

  ggplot(trace_df, aes(x=k, y=trace, col=method, group=method)) + geom_line() +
    coord_cartesian(ylim = c(0, quantile(trace_df$trace, 0.995))) +
    scale_x_continuous(name='Step') +
    scale_y_continuous(name='Trace') +
    scale_color_discrete(name='Method')
}

#' Plot ESJD over Proposed Potential Energy
#'
#' @param hmc_res A hmc result object.
#' @export
plot_esjd_potential <- function(hmc_res) {
  plot(hmc_res$energy[-1,3], get_esjd(hmc_res)$esjd,
       col=rgb(((101:500)-101)/(500-101),0.3, 0.3), pch=16,
       xlab='Potential Energy at Sample Proposal', ylab='ESJD')
}


#' Plot ESJD over Sample Covariance Matrix Trace
#'
#' @param hmc_lst A names list of hmc result objects.
#' @export
plot_esjd_trace <- function(hmc_lst) {
  traces <- lapply(hmc_lst, function(hmc_res) {
    if (length(hmc_res$eig_values)==1) { #If length=1, eig_values contains a single NA entry.
      trace_increment <- unlist(lapply(hmc_res$mass_matrices, function(M) sum(diag(solve(M)))))
      data.frame(k = 2:length(trace_increment), trace=trace_increment[-1], esjd=get_esjd(hmc_res)$esjd)
    } else {
      reg_eig_df <- get_reg_eig(hmc_res, 'minmax')
      reg_eig_df <- aggregate(reg_eig_df$ev, by=list(reg_eig_df$k), sum)
      colnames(reg_eig_df) <- c('k', 'trace')
      tmp_esjd <- get_esjd(hmc_res)
      reg_eig_df$esjd <- tmp_esjd[(nrow(tmp_esjd)-nrow(reg_eig_df)+1):nrow(tmp_esjd), 'esjd']
      reg_eig_df
    }
  })

  trace_df <- do.call(rbind, traces)
  rownames(trace_df) <- NULL

  trace_df$method <- rep(names(hmc_lst), unlist(lapply(traces, nrow)))

  ggplot(trace_df, aes(x=trace, y=esjd, col=k)) + geom_point() +
           coord_cartesian(xlim = c(quantile(trace_df$trace, 0.001),
                                    quantile(trace_df$trace, 0.995))) +
    facet_grid(trace_df$method~.) +
    scale_x_continuous(name='Step') +
    scale_y_continuous(name='Trace')
}


