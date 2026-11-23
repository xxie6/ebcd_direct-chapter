cov_fit <- function(covmat, ebnm_fn = ebnm::ebnm_point_laplace, Kmax = 1000, verbose.lvl = 0, backfit = TRUE) {
  fl <- flash_init(covmat, var_type = 0) %>%
    flash_set_verbose(verbose.lvl) %>%
    flash_greedy(ebnm_fn = ebnm_fn, Kmax = Kmax)
  if (backfit == TRUE){
    fl <- flash_backfit(fl)
  }
  s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
  s2_diff <- Inf
  while(s2 > 0 && abs(s2_diff - 1) > 1e-4) {
    covmat_minuss2 <- covmat - diag(rep(s2, ncol(covmat)))
    fl <- fl |>
      flash_update_data(covmat_minuss2) |>
      flash_set_verbose(verbose.lvl) |>
      flash_backfit() |>
      flash_nullcheck()
    old_s2 <- s2
    s2 <- max(0, mean(diag(covmat) - diag(fitted(fl))))
    s2_diff <- s2 / old_s2
    print(s2_diff)
  }
  
  return(list(fl=fl, s2 = s2))
}
