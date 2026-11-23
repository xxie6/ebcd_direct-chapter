library(lpSolve)

compute_crossprod_similarity2 <- function(est, truth){
  K_est <- ncol(est)
  K_truth <- ncol(truth)
  n <- nrow(est)
  
  #if estimates don't have same number of columns, try padding the estimate with zeros and make cosine similarity zero
  if (K_est < K_truth){
    est <- cbind(est, matrix(rep(0, n*(K_truth-K_est)), nrow = n))
  }
  
  if (K_est > K_truth){
    truth <- cbind(truth, matrix(rep(0, n*(K_est - K_truth)), nrow = n))
  }
  
  #normalize est and truth
  norms_est <- apply(est, 2, function(x){sqrt(sum(x^2))})
  norms_est[norms_est == 0] <- Inf
  
  norms_truth <- apply(truth, 2, function(x){sqrt(sum(x^2))})
  norms_truth[norms_truth == 0] <- Inf
  
  est_normalized <- t(t(est)/norms_est)
  truth_normalized <- t(t(truth)/norms_truth)
  
  #compute matrix of cosine similarities
  cosine_sim_matrix <- abs(crossprod(est_normalized, truth_normalized))
  
  assignment_problem <- lp.assign(t(cosine_sim_matrix), direction = "max")
  return((1/K_truth)*assignment_problem$objval)
}

result <- compute_crossprod_similarity2(est, truth)