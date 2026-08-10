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

An object of class `"mnTS"` with components including:

- B0:

  Estimated intercepts (1 x n matrix).

- C:

  Estimated interaction matrix (n x n).

- B:

  Estimated covariate coefficients (p x n matrix), or `NULL`.

- sigma:

  Estimated process noise scale.

- V:

  Estimated process noise covariance matrix.

- dispersion:

  Estimated observation dispersion parameter.

- logLik:

  Log-likelihood at the maximum.

- AIC:

  Akaike Information Criterion.

- mu:

  Fitted multinomial probabilities (Tmax x n matrix).

- y:

  Fitted latent states on the log-ratio scale (Tmax x n matrix).

- se.y.fitted, se.mu.upper.fitted, se.mu.lower.fitted:

  Kalman-filter standard errors and probability bands.

- par:

  Estimated parameters as a named vector.

- se:

  Approximate standard errors from the observed information matrix.

- inv.information.matrix:

  Inverse of the observed information matrix, or `NA` if not invertible.

## See also

[`mnGLMM()`](https://quinnasena.github.io/multinomialTS/reference/mnGLMM.md)
for obtaining starting values,
[`bootstrap.mnTS()`](https://quinnasena.github.io/multinomialTS/reference/bootstrap.mnTS.md)
for bootstrapped standard errors, and
[`simulate.mnTS()`](https://quinnasena.github.io/multinomialTS/reference/simulate.mnTS.md).

## Examples

``` r
# Tsample gives the time steps at which Y was observed: Y holds only the
# sampled rows, while X must cover every time step the model iterates over.
# This is a small, fast example (reference taxon plus two others over the
# first 60 time steps). See the package vignette for a full analysis.
data(story_pollen_matrix)
data(story_char_matrix)

Y_all <- story_pollen_matrix[, 1:3]
sample_idx <- which(rowSums(Y_all) != 0)
sample_idx <- sample_idx[sample_idx <= 60]
Y <- Y_all[sample_idx, ]
X <- scale(story_char_matrix)[1:60, , drop = FALSE]
n <- ncol(Y)

# NA = estimate, a number = hold fixed at that value.
B0.fixed <- matrix(c(0, rep(NA, n - 1)), nrow = 1, ncol = n)
B.fixed <- matrix(NA, nrow = ncol(X), ncol = n)
B.fixed[, 1] <- 0                    # reference taxon
V.fixed <- matrix(NA, n, n)
V.fixed[1] <- 1                      # identifiability constraint
C.start <- 0.5 * diag(n)             # self-regulation only
C.fixed <- C.start
C.fixed[C.fixed != 0] <- NA

# \donttest{
ts_mod <- mnTS(Y = Y, X = X, Tsample = sample_idx,
               B0.fixed = B0.fixed, B.fixed = B.fixed,
               C.start = C.start, C.fixed = C.fixed,
               V.fixed = V.fixed, dispersion.fixed = 1)
summary(ts_mod)
#> 
#> Call: mnTS with Tmax =  32  n =  3 
#> 
#> logLik = 93.27879,  AIC = -160.5576 [df = 13]
#> 
#> dispersion parameter =  1
#> 
#> Fitted Coefficients with approximate se
#>                                               Coef.         se           t
#> hardwood                                0.359905639 0.07141785  5.03943503
#> Fagus grandifolia                       0.197298746 0.65764588  0.30000758
#> sp.other.other                          0.908636931 0.05897038 15.40836080
#> sp.hardwood.hardwood                    0.889980048 0.09023657  9.86274249
#> sp.Fagus grandifolia.Fagus grandifolia  0.956355954 0.04445309 21.51382358
#> char_acc.hardwood                      -0.005925982 0.06924150 -0.08558425
#> char_acc.Fagus grandifolia             -0.177209792 0.12367134 -1.43290914
#>                                                    P
#> hardwood                                4.669081e-07
#> Fagus grandifolia                       7.641714e-01
#> sp.other.other                          1.438137e-53
#> sp.hardwood.hardwood                    6.037731e-23
#> sp.Fagus grandifolia.Fagus grandifolia 1.155740e-102
#> char_acc.hardwood                       9.317969e-01
#> char_acc.Fagus grandifolia              1.518838e-01
#> 
#> 
#> Overall model
#> 
#> B0 =  
#>             other  hardwood Fagus grandifolia
#> (intercept)     0 0.3599056         0.1972987
#> 
#> 
#> B =  
#>          other     hardwood Fagus grandifolia
#> char_acc     0 -0.005925982        -0.1772098
#> 
#> 
#> C =  
#>                       other hardwood Fagus grandifolia
#> other             0.9086369  0.00000          0.000000
#> hardwood          0.0000000  0.88998          0.000000
#> Fagus grandifolia 0.0000000  0.00000          0.956356
#> 
#> sigma =  0.346374
#> 
#> V =  
#>                       other  hardwood Fagus grandifolia
#> other             1.0000000 0.5139845         0.5463239
#> hardwood          0.5139845 0.3014088         0.1457267
#> Fagus grandifolia 0.5463239 0.1457267         0.9588457
# }
```
