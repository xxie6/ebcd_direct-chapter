library(RSpectra)
pca.wrapper <- function(input, args){
  Kmax <- ceiling(args$Kmax_factor * input$K)
  S_pca <- eigs_sym(input$YYt, Kmax)
  est_L <- S_pca$vectors %*% diag(sqrt(S_pca$values))
  return(list(pca_fit = S_pca, est_L = est_L, est_LLt = tcrossprod(est_L)))
}
pca_data <- pca.wrapper(input, args)
