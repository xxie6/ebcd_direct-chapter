library(flashier)
library(dplyr)
source("ebmfcov_diag_v2.R")
ebmfcov_diag.wrapper <- function(input, args){
  Kmax <- ceiling(args$Kmax_factor * input$K)
  ebmfcov_diag_fit <- cov_fit(covmat = input$YYt, 
                             ebnm_fn = args$ebnm_fn,
                             Kmax = Kmax,
                             backfit = args$backfit)
  flash_cov_ldf <- ldf(ebmfcov_diag_fit$fl)
  flash_cov_L <- flash_cov_ldf$L %*% diag(sqrt(flash_cov_ldf$D))
  return(list(ebmfcov_diag_fit = ebmfcov_diag_fit, 
              est_LLt = tcrossprod(flash_cov_L)))
}
ebmfcov_diag_data <- ebmfcov_diag.wrapper(input, args)
