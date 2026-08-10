# Wrappers around the compiled likelihoods -------------------------------------
#
# The *_ml_cpp_list functions themselves are generated into R/RcppExports.R by
# Rcpp::compileAttributes(); these wrappers add the R-side argument handling and
# expose the scalar log-likelihood the optimisers need.

#' @useDynLib multinomialTS, .registration = TRUE
#' @import Rcpp
#' @noRd
mnTS_ml_cpp_listout <- function(par, par_fixed, Y, X, Tsample) {
  # RcppArmadillo cannot coerce R's NULL to an arma::mat, so a model with no
  # covariates would fail at this boundary with "Not a matrix." before reaching
  # the no-covariate branches. An empty matrix is how the C++ side signals
  # "no covariates" (see X.is_empty() in src/source_mnTS_list.cpp).
  if (is.null(X)) X <- matrix(0, 0, 0)
  mnTS_ml_cpp_list(par, par_fixed, Y, X, Tsample)
}

#' @noRd
mnTS.ml.wrapper <- function(par, par.fixed, Y, X, Tsample) {
  mnTS_ml_cpp_listout(par, par.fixed, Y, X, Tsample)$LL
}


#' @noRd
mnGLMM_ml_cpp_listout <- function(par, par_fixed, Y, X) {
  mnGLMM_ml_cpp_list(par, par_fixed, Y, X)
}

#' @noRd
mnGLMM.ml.wrapper <- function(par, par.fixed, Y, X) {
  mnGLMM_ml_cpp_listout(par, par.fixed, Y, X)$LL
}
