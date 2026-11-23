library(flashier)
library(dplyr)
source("laplace_split_ebmfcov_diag.R")
ls_ebmfcov_diag.wrapper <- function(input, args){
  Kmax <- ceiling(args$Kmax_factor * input$K)
  ls_ebmfcov_diag_fit <- laplace_split_cov_fit(covmat = input$YYt, 
                             ebnm_fn = args$ebnm_fn,
                             Kmax = Kmax)
  flash_cov_ldf <- ldf(ls_ebmfcov_diag_fit$fl)
  flash_cov_L <- flash_cov_ldf$L %*% diag(sqrt(flash_cov_ldf$D))
  return(list(ls_ebmfcov_diag_fit = ls_ebmfcov_diag_fit, 
              est_LLt = tcrossprod(flash_cov_L)))
}
ls_ebmfcov_diag_data <- ls_ebmfcov_diag.wrapper(input, args)