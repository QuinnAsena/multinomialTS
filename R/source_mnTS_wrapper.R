# mnTS Wrapper ----------------------------------------------------------------
#' mnTS Wrapper
#'
#' This is a simple function that, by default, prints "Hello world". You can
#' customize the text to print
#'
#' @param par something
#' @param par.fixed More things
#' @param Y inout matrix of
#' @param X input matrix of
#'
#' @return This function returns a phrase to print, with or without an
#'    exclamation point added. As a side effect, this function also prints out
#'    the phrase.
#'
#'
#' @export

mnTS.ml.wrapper <- function(par, par.fixed, Y, X, Tsample) {
  retval <- mnTS_ml_cpp_listout(par, par.fixed, Y, X, Tsample)
  return(retval$LL)
}


#' mnTS_ml_cpp_list
#'
#' This is a description of what the function does.
#' @useDynLib multinomialTS, .registration = TRUE
#' @import Rcpp
#' @import RcppArmadillo
#' @export
mnTS_ml_cpp_listout <- function(par, par_fixed, Y, X, Tsample) {
  .Call('_multinomialTS_mnTS_ml_cpp_list', par, par_fixed, Y, X, Tsample, PACKAGE = 'multinomialTS')
}
