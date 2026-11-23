laplace_split_initialization <- function(S, Kmax, verbose = 2, backfit_maxiter = 500, backfit_tol = NULL){
  # change backfit_tol = NULL as default setting
  # backfitting is important in the tree setting
  
  # fit point-laplace fit with flash
  flash_laplace_fit <- flash_init(data = S, var_type = 0) |>
    flash_set_verbose(verbose = verbose) |>
    flash_greedy(Kmax = Kmax, ebnm_fn = ebnm::ebnm_point_laplace) |>
    flash_backfit(maxiter = backfit_maxiter, tol = backfit_tol) |>
    flash_nullcheck()
  
  # rescale fit so that L and F are of the same scale
  flash_laplace_fit_scaled <- ldf(flash_laplace_fit, 'i')
  LL <- flash_laplace_fit_scaled$L
  ##FF <- flash_laplace_fit_scaled$F
  
  # split into positive and negative components
  LL <- cbind(pmax(LL, 0), pmax(-LL, 0))
  ##FF <- cbind(pmax(FF, 0), pmax(-FF, 0))
  
  # remove columns of zeros
  idx.nonzero <- apply(LL, 2, function(x){return(sum(x^2))}) > 10^(-10)
  LL <- LL[, idx.nonzero]
  
  # refit weights by least squares
  n <- nrow(S)
  llt_vec <- matrix(rep(0, ncol(LL)*n*n), ncol = ncol(LL))
  for (i in 1:ncol(LL)){
    llt_vec[,i] <- c(LL[,i]%*%t(LL[,i])) 
  }
  
  nnlm_fit <- NNLM::nnlm(llt_vec, as.matrix(c(S), ncol = 1), alpha = c(0,0,0))
  
  indices_keep <- (nnlm_fit$coefficients > 0)
  LL_scaled <- LL[,indices_keep] %*% diag(sqrt(nnlm_fit$coefficients[indices_keep]))
  return(LL_scaled)
}

laplace_split_cov_fit <- function(covmat, ebnm_fn = ebnm::ebnm_point_laplace, Kmax = 1000, 
                    verbose.lvl = 0, backfit_maxiter = 500, 
                    backfit_tol = 3.81e-04) {
  # laplace-split initialization
  init_L <- laplace_split_initialization(covmat, Kmax = Kmax)

  # ebmfcov_diag fit
  fl <- flash_init(covmat, var_type = 0) |>
    flash_set_verbose(verbose.lvl) |>
    flash_factors_init(list(init_L, init_L), ebnm_fn = ebnm_fn)
  s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
  s2_diff <- Inf
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    covmat_minuss2 <- covmat - diag(rep(s2, ncol(covmat)))
    fl <- fl |>
      flash_update_data(covmat_minuss2) |>
      flash_set_verbose(verbose.lvl) |>
      flash_backfit(maxiter = backfit_maxiter, tol = backfit_tol) |>
      flash_nullcheck()
    old_s2 <- s2
    s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
  }
  
  return(list(fl=fl, s2 = s2))
}
