source("ebcd_functions_dsc.R")
ebcd.wrapper <- function(input, args){
  set.seed(args$seed)
  Kmax <- ceiling(args$Kmax_factor * input$K)
  ebcd_fit <- ebcd_nonneg(X = t(input$Y), Kmax = Kmax, maxiter_backfit = 10000, ebnm_fn = args$ebnm_fn)
  return(list(ebcd_fit = ebcd_fit, est_LLt = tcrossprod(ebcd_fit$EL)))
}
ebcd_data <- ebcd.wrapper(input, args)