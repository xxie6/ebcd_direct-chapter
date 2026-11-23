library(flashier)
library(dplyr)
flash_seminmf.wrapper <- function(input, args){
  Kmax <- ceiling(args$Kmax_factor * input$K)
  flash_fit <- flash_init(data = input$Y, var_type = 0) %>%
    flash_greedy(ebnm_fn = c(args$ebnm_fn, ebnm::ebnm_point_laplace), Kmax = Kmax)
  if (args$backfit == TRUE){
    flash_fit <- flash_backfit(flash_fit)
  }
  flash_fit <- flash_nullcheck(flash_fit)
  if (flash_fit$n_factors == 0){
    flash_fit$L_pm <- matrix(0, ncol = 1, nrow = nrow(input$Y))
    return(list(flash_fit = flash_fit, est_LLt = matrix(0, ncol = nrow(input$Y), nrow = nrow(input$Y))))
  }
  flash_ldf <- ldf(flash_fit, type = "f") # I think the default option is type = f
  flash_scaled_L <- (1/sqrt(ncol(input$Y)))*flash_ldf$L %*% diag(flash_ldf$D)
  return(list(flash_fit = flash_fit, est_LLt = tcrossprod(flash_scaled_L)))
}
flash_seminmf_data <- flash_seminmf.wrapper(input, args)
