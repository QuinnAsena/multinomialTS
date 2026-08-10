# Summary Method for Multinomial State-Space Models

Provides a summary of the fitted multinomial state-space model.

## Usage

``` r
# S3 method for class 'mnTS'
summary(object, ...)
```

## Arguments

- object:

  An object of class `"mnTS"`, typically returned by the
  [`mnTS`](https://quinnasena.github.io/multinomialTS/reference/mnTS.md)
  function.

- ...:

  Additional arguments (currently unused).

## Value

`object`, invisibly. Called for the side effect of printing the model
summary to the console.

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
summary(ts_mod)
#> 
#> Call: mnTS with Tmax =  32  n =  3 
#> 
#> logLik = 92.3436,  AIC = -158.6872 [df = 13]
#> 
#> dispersion parameter =  1
#> 
#> Fitted Coefficients with approximate se
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
#> 
#> 
#> Overall model
#> 
#> B0 =  
#>             other  hardwood Fagus grandifolia
#> (intercept)     0 0.3659118        0.07421799
#> 
#> 
#> B =  
#>          other    hardwood Fagus grandifolia
#> char_acc     0 -0.01060566        -0.2087406
#> 
#> 
#> C =  
#>                       other  hardwood Fagus grandifolia
#> other             0.9145758 0.0000000         0.0000000
#> hardwood          0.0000000 0.8769075         0.0000000
#> Fagus grandifolia 0.0000000 0.0000000         0.9511184
#> 
#> sigma =  0.2901503
#> 
#> V =  
#>                       other  hardwood Fagus grandifolia
#> other             1.0000000 0.6147506         0.5931027
#> hardwood          0.6147506 0.4304757         0.1648986
#> Fagus grandifolia 0.5931027 0.1648986         1.3522053
# }
```
