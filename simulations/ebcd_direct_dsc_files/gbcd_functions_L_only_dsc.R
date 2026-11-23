# Initialize the EB-NMF fit to covariance matrix YY' s.t. E[YY'] =
# LL'+ D from an estimate of L without nonnegative constraints.
init_cov_ebnmf <- function(fl, kset = 1:ncol(fl$flash_fit$EF[[1]])) {
  LL <- fl$flash_fit$EF[[1]][, kset, drop = FALSE]
  FF <- fl$flash_fit$EF[[2]][, kset, drop = FALSE]
  LL <- cbind(pmax(LL, 0), pmax(-LL, 0))
  LL <- cbind(fl$flash_fit$EF[[1]][, -kset, drop = FALSE], LL)
  FF <- cbind(pmax(FF, 0), pmax(-FF, 0))
  FF <- cbind(fl$flash_fit$EF[[2]][, -kset, drop = FALSE], FF)
  to.keep <- (colSums(LL) > .Machine$double.eps) &
    (colSums(FF) > .Machine$double.eps)
  LL <- LL[, to.keep, drop = FALSE]
  FF <- FF[, to.keep, drop = FALSE]
  return(list(LL, FF))
}


# Fit EBMF to covariance matrix YY' s.t. E[YY'] = LL'+ D, where D =
# s2*I and I is an identity matrix.
#
fit_ebmf_to_YY <- function (dat, fl, extrapolate = TRUE, warmstart = TRUE,
                            maxiter = 500, tol = NULL, epsilon = 2e-2,
                            verbose = 1) {
  if (is.matrix(dat) || inherits(dat, "Matrix")) {
    data_diag <- diag(dat)
  } else {
    data_diag <- Matrix::rowSums(dat$U * dat$D * dat$V)
  }
  
  s2 <- max(0, mean(data_diag - rowSums(fl$L_pm * fl$F_pm)))
  s2_diff <- Inf
  old_s2 <- 0
  
  ### alternate between estimating s2 and backfitting until convergence
  while(s2 > 0 && abs(s2_diff - 1) > epsilon) {
    if (is.matrix(dat) || inherits(dat, "Matrix")) {
      dat_minuss2 <- dat - diag(rep(s2, ncol(dat)))
    }
    else{
      dat_minuss2 <- list(U = dat$U, D = dat$D, V = dat$V,
                          S = Matrix::Diagonal(nrow(dat$U), -s2))
    }
    Y2_diff <- sum((data_diag - s2)^2 - (data_diag - old_s2)^2)
    fl <- fl |>
      flash_update_data(dat_minuss2, Y2_diff = Y2_diff) |>
      flash_backfit(extrapolate = extrapolate, warmstart = warmstart,
                    maxiter = maxiter, tol = tol, verbose = verbose)
    
    old_s2 <- s2
    s2 <- max(0, mean(data_diag - rowSums(fl$L_pm * fl$F_pm)))
    s2_diff <- s2 / old_s2
  }
  
  return(list(dat = dat, fl = fl, s2 = s2))
}

fit_gbcd_L_only <- function (Y, Kmax, prior = ebnm::ebnm_generalized_binary, 
                      maxiter1 = 500, maxiter2 = 200, maxiter3 = 500,
                      control = list(), verbose = 1) {
  
  control <- modifyList(fit_gbcd_control_default(), control, keep.null = TRUE)
  extrapolate <- control$extrapolate
  warmstart <- control$warmstart
  corr_thres <- control$corr_thres
  
  ### form the covariance matrix from the cell by gene matrix of gene expression data
  print("Form cell by cell covariance matrix...")
  start_time = proc.time()
  if (2 * ncol(Y) * mean(Y > 0) < nrow(Y)) {
    # Use lowrank plus sparse representation:
    dat <- list(U = Y, D = rep(1 / ncol(Y), ncol(Y)), V = Y)
  } else {
    # Form covariance matrix:
    dat <- Matrix::tcrossprod(Y) / ncol(Y)
  }
  fit.init <- flash_init(dat, var_type = 0)
  runtime = proc.time() - start_time
  print(runtime)
  
  ### fit EBMF with point laplace prior to covariance matrix without considering the diagonal component
  print("Initialize GEP membership matrix L...")
  start_time = proc.time()
  fit.cov <- fit.init |>
    flash_greedy(Kmax = 1, ebnm_fn = ebnm_point_laplace) |>
    flash_greedy(Kmax = Kmax - 1, ebnm_fn = ebnm_point_laplace) |>
    flash_backfit(maxiter = 25, verbose = verbose)
  
  ### fit EBMF with point laplace prior to covariance matrix with the diagonal component
  fit.cov <- fit_ebmf_to_YY(dat = dat, fl = fit.cov, maxiter = maxiter1, verbose = verbose)$fl
  runtime = proc.time() - start_time
  print(runtime)
  
  ### initialize EB-NMF fit from the EBMF fit with point laplace prior
  print("Estimate GEP membership matrix L...")
  start_time = proc.time()
  cov.init <- init_cov_ebnmf(fit.cov)
  kmax <- which.max(colSums(cov.init[[1]]))
  fit.cov <- fit.init |>
    flash_factors_init(
      init = lapply(cov.init, function(x) x[, kmax, drop = FALSE]),
      ebnm_fn = ebnm_point_laplace
    ) |>
    flash_factors_init(
      init = lapply(cov.init, function(x) x[, -c(kmax), drop = FALSE]),
      ebnm_fn = prior
    ) |>
    flash_backfit(extrapolate = FALSE, warmstart = warmstart, maxiter = 25, verbose = verbose)
  
  ### keep at most Kmax factors based on proportion of variance explained and refit EB-NMF to covariance matrix
  kset <- (fit.cov$pve > 0)
  kall <- 1:fit.cov$n_factors
  if(!all(kset))
    fit.cov <- flash_factors_remove(fit.cov, kset=kall[!kset])
  fit.cov <- fit_ebmf_to_YY(dat = dat, fl = fit.cov, extrapolate = extrapolate, warmstart = warmstart, maxiter = maxiter2, verbose = verbose)$fl
  runtime = proc.time() - start_time
  print(runtime)
  
  ### I had a case where I had an additional zero column in L after running fit_ebmf_to_YY
  ### check for factors that are zero and remove them
  fit.cov <- flash_nullcheck(fit.cov)
  
  ### compute rescaled version of L
  #fit.gbcd.rescale.ldf <- ldf(fit.cov$fl, type = 'i')
  fit.gbcd.rescale.ldf <- ldf(fit.cov, type = 'i')
  fit.gbcd.rescale.L <- fit.gbcd.rescale.ldf$L %*% diag(sqrt(fit.gbcd.rescale.ldf$D))
  
  ### estimate GEP signatures by fitting EB-SNMF to gene expression data matrix with fixed L estimated from covariance decomposition above
  #print("Estimate GEP signature matrix F...")
  #start_time = proc.time()
  #res0 <- fit_ebmf_to_Y(Y, fit.cov, corr_thres, maxiter3)
  #res <- res0$gbcd_output
  #runtime = proc.time() - start_time
  #print(runtime)
  
  k.order <- order(fit.cov$pve, decreasing = TRUE)
  fit.L <- fit.cov$L_pm[, k.order]
  fit.L <- t(t(fit.L)/apply(fit.L, 2, max))
  corr <- diag(cor(fit.cov$L_pm[, k.order], fit.cov$F_pm[, k.order]))
  k.idx <- corr > corr_thres
  L.pm <- fit.L[, k.idx]
  colnames(L.pm) <- c("Baseline", paste0("GEP", 1:(ncol(L.pm)-1)))
  res <- list(L = L.pm)
  
  
  ### only keep columns of rescaled L with high correlation
  fit.gbcd.rescale.L <- fit.gbcd.rescale.L[,k.order][,k.idx]
  
  return(list(res = res, scaled_L = fit.gbcd.rescale.L))
}

fit_gbcd_control_default <- function()
  list(warmstart = TRUE,
       extrapolate = FALSE,
       corr_thres = 0.8)