source("sindclus.R")
sindclus.wrapper <- function(input, args){
  Kmax <- ceiling(args$Kmax_factor * input$K)
  set.seed(args$seed)
  sindclus_fit <- fit_sindclus(input$YYt, 
                               Kmax, 
                               args$init_method, 
                               conv_tol = 1e-6, 
                               args$additive_term, 
                               args$off_diagonal)
  return(list(sindclus_fit = sindclus_fit,
              est_LLt = (sindclus_fit$P %*% diag(sindclus_fit$W) %*% t(sindclus_fit$P))))
}
sindclus_data <- sindclus.wrapper(input, args)
