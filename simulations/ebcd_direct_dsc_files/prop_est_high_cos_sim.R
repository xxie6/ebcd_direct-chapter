compute_prop_est_high_cos_sim <- function(est, truth, args){
  K_est <- ncol(est)
  
  norms_est <- apply(est, 2, function(x){sqrt(sum(x^2))})
  norms_est[norms_est == 0] <- Inf
  
  norms_truth <- apply(truth, 2, function(x){sqrt(sum(x^2))})
  norms_truth[norms_truth == 0] <- Inf
  
  est_normalized <- t(t(est)/norms_est)
  truth_normalized <- t(t(truth)/norms_truth)
  
  #compute matrix of cosine similarities
  cosine_sim_matrix <- abs(crossprod(est_normalized, truth_normalized))
  
  est_highest_cos_sim <- apply(cosine_sim_matrix, 1, max)
  num_est_high_cos_sim <- sum(est_highest_cos_sim > args$threshold)
  return(num_est_high_cos_sim/K_est)
}

result <- compute_prop_est_high_cos_sim(est, truth, args)