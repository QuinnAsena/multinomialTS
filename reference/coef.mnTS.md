# Coefficients for Multinomial State-Space Model

Extracts estimated coefficients from a fitted multinomial state-space
model.

## Usage

``` r
# S3 method for class 'mnTS'
coef(object, ...)
```

## Arguments

- object:

  An object of class `"mnTS"`, as returned by the
  [`mnTS`](https://quinnasena.github.io/multinomialTS/reference/mnTS.md)
  function.

- ...:

  Additional arguments (currently unused).

## Value

A matrix of estimated coefficients. When an information matrix is
available it has columns `Coef.`, `se`, `t` and `P`; otherwise it holds
point estimates only.

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
coef(ts_mod)
#>                                              Coef.         se           t
#> hardwood                                0.36591181        NaN         NaN
#> Fagus grandifolia                       0.07421799 0.64074684  0.11583044
#> sp.other.other                          0.91457579        NaN         NaN
#> sp.hardwood.hardwood                    0.87690755        NaN         NaN
#> sp.Fagus grandifolia.Fagus grandifolia  0.95111844 0.02747347 34.61952616
#> char_acc.hardwood                      -0.01060566 0.11612721 -0.09132796
#> char_acc.Fagus grandifolia             -0.20874064 0.15047134 -1.38724516
#>                                                    P
#> hardwood                                         NaN
#> Fagus grandifolia                       9.077869e-01
#> sp.other.other                                   NaN
#> sp.hardwood.hardwood                             NaN
#> sp.Fagus grandifolia.Fagus grandifolia 1.284612e-262
#> char_acc.hardwood                       9.272320e-01
#> char_acc.Fagus grandifolia              1.653670e-01
# }
```
