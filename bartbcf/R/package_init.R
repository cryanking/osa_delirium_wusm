#' @useDynLib bartbcf
#' @importFrom Rcpp sourceCpp evalCpp
#' @import BART
#' @importFrom stats ar lm plogis pnorm qchisq qlogis residuals sd
NULL



#' Compute causal effect estimates using Bayesian Adaptive Regression Trees. The methods here are those of Hahn in (citation). BART tree-based estimates under treatment and control are computed and returned for each subject. Because of an implied first "split" on the treatment variable estimates of treatment effect are allowed to be highly heterogenious and ITEs and subpopulation effects are directly computable.
"_PACKAGE"