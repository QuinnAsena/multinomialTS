# Print Method for Multinomial GLMMs

Prints a fitted multinomial GLMM, using the same layout as
[`summary.mnGLMM`](https://quinnasena.github.io/multinomialTS/reference/summary.mnGLMM.md).

## Usage

``` r
# S3 method for class 'mnGLMM'
print(x, ...)
```

## Arguments

- x:

  An object of class `"mnGLMM"`, as returned by the
  [`mnGLMM`](https://quinnasena.github.io/multinomialTS/reference/mnGLMM.md)
  function.

- ...:

  Additional arguments passed to
  [`summary.mnGLMM`](https://quinnasena.github.io/multinomialTS/reference/summary.mnGLMM.md).

## Value

`x`, invisibly.

## Examples

``` r
data(story_pollen_matrix)
sample_idx <- which(rowSums(story_pollen_matrix) != 0)
n <- ncol(story_pollen_matrix)
glmm_mod <- mnGLMM(Y = story_pollen_matrix[sample_idx, ], V.fixed = diag(n))
print(glmm_mod)
#> 
#> Call: mnGLMM with Tmax =  93  n =  5 
#> 
#> logLik = 449.3075,  AIC = -888.6149 [df = 5]
#> 
#> dispersion parameter =  1
#> 
#> Fitted Coefficients with approximate se
#>                                     Coef.         se         t            P
#> (intercept).hardwood           0.09345372 0.07082625  1.319478 1.870092e-01
#> (intercept).Fagus grandifolia -0.31836611 0.07195283 -4.424650 9.659870e-06
#> (intercept).Ulmus             -0.39757671 0.07161670 -5.551453 2.833053e-08
#> (intercept).Quercus            0.84814202 0.07232762 11.726392 9.335386e-32
#> 
#> 
#> Overall model
#> 
#> B =  
#>      [,1]       [,2]       [,3]       [,4]     [,5]
#> [1,]    0 0.09345372 -0.3183661 -0.3975767 0.848142
#> 
#> sigma =  0.495335
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
