compute_diff_num_factors <- function(est, truth){
  K_est <- ncol(est)
  K_truth <- ncol(truth)
  if ((K_est == 1) & (sum(K_est) == 0)){
    K_est <- 0
  }
  K_diff <- K_est - K_truth
  return(K_diff)
}

result <- compute_diff_num_factors(est, truth)
