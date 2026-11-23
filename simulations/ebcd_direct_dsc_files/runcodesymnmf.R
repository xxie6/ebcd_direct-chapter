library(codesymnmf)
codesymnmf.wrapper <- function(input, args){
  Kmax <- ceiling(args$Kmax_factor * input$K)
  codesymnmf_results <- codesymnmf(input$YYt, r = Kmax)
  return(list(codesymnmf_fit = codesymnmf_results, est_LLt = tcrossprod(codesymnmf_results$H)))
}
codesymnmf_data <- codesymnmf.wrapper(input, args)