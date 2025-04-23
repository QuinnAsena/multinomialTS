# mnTS Wrapper ----------------------------------------------------------------
#' mnTS Wrapper
#'
#' These are internal functions required by multinomialTS.

mnTS.ml.wrapper <- function(par, par.fixed, Y, X, Tsample) {
  retval <- mnTS_ml_cpp_listout(par, par.fixed, Y, X, Tsample)
  return(retval$LL)
}


#' mnTS_ml_cpp_list
#'
#' These are internal functions required by multinomialTS.
#' @useDynLib multinomialTS, .registration = TRUE
#' @import Rcpp
#' @import RcppArmadillo

mnTS_ml_cpp_listout <- function(par, par_fixed, Y, X, Tsample) {
  .Call('_multinomialTS_mnTS_ml_cpp_list', par, par_fixed, Y, X, Tsample, PACKAGE = 'multinomialTS')
}
