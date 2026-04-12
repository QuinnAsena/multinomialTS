# Multinomial Generalized Linear Mixed Model (mnGLMM)

Fits a multivariate generalized linear regression model for
multinomially distributed response data. This function estimates
driver-species relationships (B coefficients), observation-level
dispersion, and species covariance structures via maximum likelihood or
REML.

## Usage

``` r
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
```

## Arguments

- Y:

  A matrix of multinomially distributed count data (e.g., community
  count data).

- X:

  A matrix of covariates (predictors), which may be of mixed type
  Covariates should be scaled when appropriate. Can be `NULL`.

- B.fixed:

  A matrix indicating which B coefficients (driver-species
  relationships) to estimate. The number of columns must equal `ncol(Y)`
  (number of species), and the number of rows must equal `ncol(X) + 1`
  (intercept + covariates).

- B.start:

  A matrix of starting values for the B coefficients. Dimensions must
  match `B.fixed`.

- sigma.fixed:

  Fixed value for the overall model variance. Use `NA` to estimate it
  from the model.

- sigma.start:

  Starting value for estimating `sigma.fixed` (default is 0.1).

- dispersion.fixed:

  Fixed dispersion parameter to account for over- or under-dispersion. A
  value of 1 corresponds to no extra dispersion (pure multinomial).

- dispersion.start:

  Starting value for estimating `dispersion.fixed`.

- V.fixed:

  A species-by-species covariance matrix representing environmental
  variation.

- V.start:

  Starting values for `V.fixed`.

- method:

  Optimization method. Acceptable values include `"Nelder-Mead"`,
  `"BFGS"` (via `optim`), and `"bobyqa"` (via the `minqa` package).

- optim.control:

  Optional list of control parameters passed to the optimizer. See the
  `minqa` package documentation for details.

- maxit.optim:

  Maximum number of iterations for the optimizer (default is 1e+05).
  Increase if of optimizer needs more itterations.

- REML:

  Logical. If `TRUE`, uses restricted maximum likelihood for variance
  estimation.

- compute.information.matrix:

  Logical. If `TRUE`, computes the observed information matrix.

## Value

An object of class `"mnGLMM"` with components including:

- B:

  Estimated covariate coefficients (p x n matrix).

- sigma:

  Estimated process noise scale.

- V:

  Estimated species covariance matrix.

- dispersion:

  Estimated observation dispersion parameter.

- mu:

  Fitted multinomial probabilities (Tmax x n matrix).

- logLik:

  Log-likelihood at the maximum.

- AIC:

  Akaike Information Criterion.

- par:

  Estimated parameters as a named vector.

- se:

  Approximate standard errors from the observed information matrix.

- inv.information.matrix:

  Inverse of the observed information matrix, or `NA` if not invertible.
