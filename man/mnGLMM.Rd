% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/source_mnGLMM.R
\name{mnGLMM}
\alias{mnGLMM}
\title{Multinomial GLMM}
\usage{
mnGLMM(
  Y,
  X = NULL,
  B.fixed = if (is.null(X)) matrix(c(0, rep(NA, ncol(Y) - 1)), nrow = 1, ncol = ncol(Y))
    else matrix(c(rep(0, (ncol(X) + 1) * ncol(Y)), rep(NA, (ncol(Y) - 1) * (ncol(X) +
    1))), nrow = ncol(X) + 1, ncol = ncol(Y)),
  B.start = if (is.null(X)) matrix(0, nrow = 1, ncol = ncol(Y)) else matrix(0, nrow =
    ncol(X) + 1, ncol = ncol(Y)),
  sigma.fixed = NA,
  sigma.start = 0.1,
  dispersion.fixed = 1,
  dispersion.start = 1,
  V.fixed = diag(ncol(Y)),
  V.start = diag(ncol(Y)),
  method = "bobyqa",
  optim.control = NULL,
  maxit.optim = 1e+05,
  REML = FALSE,
  compute.information.matrix = TRUE
)
}
\arguments{
\item{Y}{inout matrix of}

\item{X}{input matrix of}

\item{REML}{Description}

\item{fitted.values}{TRUE/FALSE argument indicating}
}
\value{
This function returns a phrase to print, with or without an
   exclamation point added. As a side effect, this function also prints out
   the phrase.
}
\description{
This is a not a simple function
customize the text to print.
}
