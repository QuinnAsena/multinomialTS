# Bootstrap method for mnTS objects

Provides bootstrapped confidence intervals for an object of class
`mnTS`.

## Usage

``` r
# S3 method for class 'mnTS'
bootstrap(object, reps, dispersion.fixed = 1, maxit.optim = 1e+05, ...)
```

## Arguments

- object:

  An object of class `mnTS`.

- reps:

  Integer. Number of bootstrap replications. More is better; the package
  vignette discusses how to choose a value.

- dispersion.fixed:

  Numeric. Dispersion parameter (default is 1).

- maxit.optim:

  Integer. Max iterations for optimization (default 100000).

- ...:

  Additional arguments (not currently used).

## Value

A list with three elements: `coef.table.boot`, a table of bootstrapped
coefficients with standard deviations, t-scores, P-values and 68%
bounds; `all_mods_pars`, one row of coefficients per successful
replicate; and `all_mods`, the fitted replicate models.

## Examples

``` r
# A small, fast example: the reference taxon plus two others over the first
# 60 time steps. See the package vignette for a full analysis.
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
boot_out <- bootstrap(ts_mod, reps = 3)
boot_out$coef.table.boot
#>                 boot_mean           se          t            P       upper68
#> y2           2.176641e-01 5.214389e-02  4.1742976 2.989069e-05  2.698080e-01
#> y3          -1.379883e-01 5.645390e-01 -0.2444266 8.069005e-01  4.265507e-01
#> sp.y1.y1     1.177357e-01 3.710052e-01  0.3173425 7.509837e-01  4.887409e-01
#> sp.y2.y2     9.497331e-01 1.705068e-02 55.7005925 0.000000e+00  9.667838e-01
#> sp.y3.y3     8.716086e-01 2.134049e-01  4.0842948 4.421087e-05  1.085014e+00
#> v.y1.y2      5.577487e-01 7.154506e-01  0.7795768 4.356400e-01  1.273199e+00
#> v.y2.y2      3.106642e-01 2.824374e-01  1.0999398 2.713584e-01  5.931016e-01
#> v.y1.y3      2.065884e-01 1.101435e+00  0.1875630 8.512192e-01  1.308023e+00
#> v.y2.y3     -4.180908e-01 5.332095e-01 -0.7841022 4.329801e-01  1.151188e-01
#> v.y3.y3     -2.383577e-06 4.094654e-06 -0.5821192 5.604864e-01  1.711077e-06
#> char_acc.y2  9.680478e-03 9.391296e-02  0.1030793 9.179001e-01  1.035934e-01
#> char_acc.y3 -3.065430e-01 7.630476e-02 -4.0173507 5.885611e-05 -2.302382e-01
#>                   lower68 n_fail_converged n_total
#> y2           1.655202e-01                0       3
#> y3          -7.025273e-01                0       3
#> sp.y1.y1    -2.532695e-01                0       3
#> sp.y2.y2     9.326825e-01                0       3
#> sp.y3.y3     6.582037e-01                0       3
#> v.y1.y2     -1.577019e-01                0       3
#> v.y2.y2      2.822673e-02                0       3
#> v.y1.y3     -8.948462e-01                0       3
#> v.y2.y3     -9.513003e-01                0       3
#> v.y3.y3     -6.478231e-06                0       3
#> char_acc.y2 -8.423248e-02                0       3
#> char_acc.y3 -3.828477e-01                0       3
# }
```
