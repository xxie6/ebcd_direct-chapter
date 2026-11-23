library(symebcovmf)
symebcovmf.wrapper <- function(input, args){
  set.seed(args$seed)
  Kmax <- ceiling(args$Kmax_factor * input$K)
  symebcovmf_fit <- sym_ebcovmf_fit(S = input$YYt, Kmax = Kmax, maxiter = 100, rank_one_tol = 10^(-5), tol = 10^(-5), ebnm_fn = args$ebnm_fn, refit_lam = TRUE)
  if (args$backfit == TRUE){
    symebcovmf_fit <- sym_ebcovmf_backfit(S = input$YYt, sym_ebcovmf_obj = symebcovmf_fit, ebnm_fn = args$ebnm_fn, 
                                          backfit_maxiter = 100, backfit_tol = 10^(-5), optim_maxiter = 100, optim_tol = 10^(-5))
  }
  est_K <- length(symebcovmf_fit$lambda)
  return(list(symebcovmf_fit = symebcovmf_fit, est_LLt = symebcovmf_fit$L_pm %*% diag(symebcovmf_fit$lambda, ncol = est_K) %*% t(symebcovmf_fit$L_pm)))
}
symebcovmf_data <- symebcovmf.wrapper(input, args)