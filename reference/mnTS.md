# Multinomial State-Space Model (mnTS)

Fits a multinomial state-space model for multivariate count data,
allowing for latent temporal processes, covariate effects, and species
interactions.

## Usage

``` r
mnTS(
  Y,
  X = NULL,
  Tsample = 1:nrow(Y),
  B0.fixed = matrix(c(0, rep(NA, ncol(Y) - 1)), nrow = 1, ncol = ncol(Y)),
  B0.start = matrix(0, nrow = 1, ncol = ncol(Y)),
  C.fixed = diag(rep(NA, ncol(Y))),
  C.start = 0.01 * diag(ncol(Y)),
  B.fixed = if (is.null(X)) NULL else matrix(NA, nrow = ncol(X), ncol = ncol(Y)),
  B.start = if (is.null(X)) NULL else matrix(0, nrow = ncol(X), ncol = ncol(Y)),
  sigma.fixed = NA,
  sigma.start = 0.1,
  dispersion.fixed = 1,
  dispersion.start = 1,
  V.fixed = diag(ncol(Y)),
  V.start = diag(ncol(Y)),
  compute.information.matrix = TRUE,
  method = "bobyqa",
  optim.control = NULL,
  maxit.optim = 1e+05,
  hessian.method.args = list(eps = 1e-04, d = 1e-04, r = 4, v = 2)
)
```

## Arguments

- Y:

  A matrix of multinomially distributed count data (e.g., community
  count data).

- X:

  A matrix of covariates (predictors), which may be of mixed type.
  Covariates should be scaled when appropriate. Can be `NULL`.

- Tsample:

  A vector of row indices specifying the subset of observations in `Y`
  to treat as temporal samples.

- B0.fixed:

  A 1 by `ncol(Y)` matrix of species intercepts to estimate.

- B0.start:

  A matrix of starting values for `B0.fixed`.

- C.fixed:

  A species-by-species matrix of interactions, indicating which
  interactions to estimate.

- C.start:

  A matrix of starting values for `C.fixed`.

- B.fixed:

  A matrix indicating which B coefficients (driver-species
  relationships) to estimate. Should have `ncol(Y)` columns and
  `ncol(X)` rows.

- B.start:

  A matrix of starting values for `B.fixed`. Dimensions of `B.start`
  should match `B.fixed`.

- sigma.fixed:

  Fixed value for the overall model variance. Use `NA` to estimate it
  from the model.

- sigma.start:

  Starting value for estimating `sigma.fixed`.

- dispersion.fixed:

  Fixed dispersion parameter for observation-level variation. A value of
  1 corresponds to no over- or under-dispersion.

- dispersion.start:

  Starting value for estimating `dispersion.fixed`.

- V.fixed:

  A species-by-species covariance matrix representing environmental
  variation.

- V.start:

  Starting values for `V.fixed`.

- compute.information.matrix:

  Logical. If `TRUE`, computes the observed information matrix.

- method:

  Optimization method. Acceptable values include `"Nelder-Mead"`,
  `"BFGS"` (via `optim`), and `"bobyqa"` (via the `minqa` package).

- optim.control:

  Optional list of control parameters passed to the optimizer. See the
  `minqa` package documentation for details.

- maxit.optim:

  Maximum number of iterations for the optimizer (default is 1e+05).
  Increase if the optimizer needs more iterations.

- hessian.method.args:

  A list of control parameters passed to the numerical Hessian
  calculator (e.g.,
  [`numDeriv::hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html)).

## Value

An object of class `"mnTS"` containing estimated parameters.
