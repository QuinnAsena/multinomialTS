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
#>                 boot_mean           se           t            P      upper68
#> y2           5.276643e-01 2.698807e-01  1.95517626 5.056227e-02 7.975449e-01
#> y3           8.668866e-01 8.050081e-01  1.07686691 2.815397e-01 1.671895e+00
#> sp.y1.y1     4.279147e-01 4.718925e-01  0.90680565 3.645096e-01 8.998072e-01
#> sp.y2.y2     8.754452e-01 1.375745e-01  6.36342590 1.973026e-10 1.013020e+00
#> sp.y3.y3     6.673375e-01 5.080871e-01  1.31343147 1.890376e-01 1.175425e+00
#> v.y1.y2      2.270793e-01 7.899540e-01  0.28745891 7.737610e-01 1.017033e+00
#> v.y2.y2      2.021561e-01 1.758526e-01  1.14957689 2.503182e-01 3.780087e-01
#> v.y1.y3      1.370357e-01 1.256096e+00  0.10909654 9.131259e-01 1.393132e+00
#> v.y2.y3     -1.470164e-02 1.139909e-01 -0.12897204 8.973798e-01 9.928930e-02
#> v.y3.y3      1.316664e-05 2.110612e-05  0.62383042 5.327390e-01 3.427276e-05
#> char_acc.y2 -3.287752e-03 2.241912e-01 -0.01466494 9.882995e-01 2.209035e-01
#> char_acc.y3  5.112835e-02 2.996864e-01  0.17060619 8.645334e-01 3.508147e-01
#>                   lower68 n_fail_converged n_total
#> y2           2.577836e-01                0       3
#> y3           6.187848e-02                0       3
#> sp.y1.y1    -4.397771e-02                0       3
#> sp.y2.y2     7.378707e-01                0       3
#> sp.y3.y3     1.592505e-01                0       3
#> v.y1.y2     -5.628747e-01                0       3
#> v.y2.y2      2.630349e-02                0       3
#> v.y1.y3     -1.119060e+00                0       3
#> v.y2.y3     -1.286926e-01                0       3
#> v.y3.y3     -7.939481e-06                0       3
#> char_acc.y2 -2.274790e-01                0       3
#> char_acc.y3 -2.485580e-01                0       3
# }
```
