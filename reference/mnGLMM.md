# Multinomial Generalized Linear Mixed Model (mnGLMM)

Fits a multivariate generalized linear regression model for
multinomially distributed response data. This function estimates
driver-species relationships (B coefficients), observation-level
dispersion, and species covariance structures via maximum likelihood.

## Usage

``` r
mnGLMM(
  Y,
  X = NULL,
  B.fixed = NULL,
  B.start = NULL,
  sigma.fixed = NA,
  sigma.start = 0.1,
  dispersion.fixed = 1,
  dispersion.start = 1,
  V.fixed = diag(ncol(Y)),
  V.start = diag(ncol(Y)),
  method = "bobyqa",
  optim.control = NULL,
  maxit.optim = 1e+05,
  compute.information.matrix = TRUE,
  hessian.method.args = list(eps = 1e-04, d = 0.1)
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
  Increase if the optimizer needs more iterations.

- compute.information.matrix:

  Logical. If `TRUE`, computes the observed information matrix.

- hessian.method.args:

  A list of control parameters passed to the numerical Hessian
  calculator (e.g.,
  [`numDeriv::hessian`](https://rdrr.io/pkg/numDeriv/man/hessian.html)).

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

## See also

[`mnTS()`](https://quinnasena.github.io/multinomialTS/reference/mnTS.md),
which uses the estimates from `mnGLMM` as starting values.

## Examples

``` r
data(story_pollen_matrix)
data(story_char_matrix)

# mnGLMM cannot handle gaps, so use only the rows with pollen counts.
sample_idx <- which(rowSums(story_pollen_matrix) != 0)
Y <- story_pollen_matrix[sample_idx, ]
X <- scale(story_char_matrix)[sample_idx, , drop = FALSE]

n <- ncol(Y)          # number of taxa
p <- ncol(X) + 1      # covariates + intercept

# Fix the reference taxon (column 1) to 0; estimate the rest (NA = estimate).
B.fixed <- matrix(c(rep(0, p), rep(NA, (n - 1) * p)), p, n)
B.start <- matrix(c(rep(0, p), rep(0.01, (n - 1) * p)), p, n)

glmm_mod <- mnGLMM(Y = Y, X = X, B.start = B.start, B.fixed = B.fixed,
                   V.fixed = diag(n))
summary(glmm_mod)
#> 
#> Call: mnGLMM with Tmax =  93  n =  5 
#> 
#> logLik = 464.0599,  AIC = -910.1198 [df = 9]
#> 
#> dispersion parameter =  1
#> 
#> Fitted Coefficients with approximate se
#>                                     Coef.         se          t            P
#> (intercept).hardwood           0.11645835 0.07016323  1.6598202 9.695062e-02
#> char_acc.hardwood             -0.12435609 0.06504119 -1.9119590 5.588146e-02
#> (intercept).Fagus grandifolia -0.32273476 0.07108291 -4.5402583 5.618534e-06
#> char_acc.Fagus grandifolia    -0.40189296 0.07939938 -5.0616638 4.156135e-07
#> (intercept).Ulmus             -0.37828175 0.07095489 -5.3312994 9.751250e-08
#> char_acc.Ulmus                -0.23302904 0.07430120 -3.1362755 1.711084e-03
#> (intercept).Quercus            0.84797318 0.07115093 11.9179493 9.542857e-33
#> char_acc.Quercus              -0.02854291 0.07145099 -0.3994754 6.895430e-01
#> 
#> 
#> Overall model
#> 
#> B =  
#>      [,1]       [,2]       [,3]       [,4]        [,5]
#> [1,]    0  0.1164584 -0.3227348 -0.3782817  0.84797318
#> [2,]    0 -0.1243561 -0.4018930 -0.2330290 -0.02854291
#> 
#> sigma =  0.4812694
#> 
#> V =  
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    0    0    0    0
#> [2,]    0    1    0    0    0
#> [3,]    0    0    1    0    0
#> [4,]    0    0    0    1    0
#> [5,]    0    0    0    0    1
#> 
```
