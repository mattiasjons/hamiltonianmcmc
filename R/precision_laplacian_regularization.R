precision_laplacian_regularization <- function(M_inv, laplace_var_explained = 0.99) {
  #M_inv <- cov(hmc_res$samples[30:40,])
  eg <- eigen(M_inv)
  lambda_shrunk <- (10 / (10 + 5.0)) * eg$values + 0.1 * (5.0 / (10 + 5.0))
  M <- eg$vectors %*% diag(1/lambda_shrunk) %*% t(eg$vectors)

  degree_matrix <- diag(rowSums(abs(M)))
  laplacian_matrix <- degree_matrix - M

  eigen_decomp <- eigen(laplacian_matrix)
  V <- eigen_decomp$vectors
  #plot(eigen_decomp$values)

  #threshold <- eigen_decomp$values[which(cumsum(eigen_decomp$values)/sum(eigen_decomp$values)>laplace_var_explained)[1]]
  n_laplacian <- 170
  threshold <- eigen_decomp$values[n_laplacian]

  lambda_sparse <- diag(ifelse(eigen_decomp$values >= threshold, eigen_decomp$values, 0))
  laplacian_matrix_sparse <- V %*% lambda_sparse %*% t(V)
  precision_matrix_sparse <- degree_matrix + laplacian_matrix_sparse

  precision_matrix_sparse
}
