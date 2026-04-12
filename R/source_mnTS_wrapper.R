# mnTS Wrapper ----------------------------------------------------------------
#' @noRd
mnTS.ml.wrapper <- function(par, par.fixed, Y, X, Tsample) {
  mnTS_ml_cpp_listout(par, par.fixed, Y, X, Tsample)$LL
}


#' @useDynLib multinomialTS, .registration = TRUE
#' @import Rcpp
#' @import RcppArmadillo
#' @noRd
mnTS_ml_cpp_listout <- function(par, par_fixed, Y, X, Tsample) {
  .Call('_multinomialTS_mnTS_ml_cpp_list', par, par_fixed, Y, X, Tsample, PACKAGE = 'multinomialTS')
}


# mnGLMM Wrappers ----------------------------------------------------------
#' @noRd
mnGLMM.ml.wrapper <- function(par, par.fixed, Y, X) {
  mnGLMM_ml_cpp_listout(par, par.fixed, Y, X)$LL
}

#' @noRd
mnGLMM_ml_cpp_listout <- function(par, par_fixed, Y, X) {
  .Call('_multinomialTS_mnGLMM_ml_cpp_list', par, par_fixed, Y, X, PACKAGE = 'multinomialTS')
}
