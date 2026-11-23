compute_prop_true_high_cos_sim <- function(est, truth, args){
  K_true <- ncol(truth)
  
  norms_est <- apply(est, 2, function(x){sqrt(sum(x^2))})
  norms_est[norms_est == 0] <- Inf
  
  norms_truth <- apply(truth, 2, function(x){sqrt(sum(x^2))})
  norms_truth[norms_truth == 0] <- Inf
  
  est_normalized <- t(t(est)/norms_est)
  truth_normalized <- t(t(truth)/norms_truth)
  
  #compute matrix of cosine similarities
  cosine_sim_matrix <- abs(crossprod(est_normalized, truth_normalized))
  
  true_highest_cos_sim <- apply(cosine_sim_matrix, 2, max)
  
  num_true_high_cos_sim <- sum(true_highest_cos_sim > args$threshold)
  return(num_true_high_cos_sim/K_true)
}

result <- compute_prop_true_high_cos_sim(est, truth, args)