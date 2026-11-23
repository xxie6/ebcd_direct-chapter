library(ebnm)
library(flashier)
library(ashr)
library(Matrix)
library(stats)
source("gbcd_functions_L_only_dsc.R")
gbcd.wrapper <- function(input,args){
  Kmax <- ceiling(args$Kmax_factor * input$K)
  gbcd_fit <- fit_gbcd_L_only(Y = input$Y, Kmax = Kmax, prior = args$ebnm_fn)
  return(list(fit_obj = gbcd_fit$res, est_LLt = tcrossprod(gbcd_fit$scaled_L)))
}
gbcd_data <- gbcd.wrapper(input,args)
