# Print Method for Multinomial State-Space Models

Prints a fitted multinomial state-space model, using the same layout as
[`summary.mnTS`](https://quinnasena.github.io/multinomialTS/reference/summary.mnTS.md).

## Usage

``` r
# S3 method for class 'mnTS'
print(x, ...)
```

## Arguments

- x:

  An object of class `"mnTS"`, as returned by the
  [`mnTS`](https://quinnasena.github.io/multinomialTS/reference/mnTS.md)
  function.

- ...:

  Additional arguments passed to
  [`summary.mnTS`](https://quinnasena.github.io/multinomialTS/reference/summary.mnTS.md).

## Value

`x`, invisibly.

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
print(ts_mod)
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
