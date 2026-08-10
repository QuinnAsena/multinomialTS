# GLMM summary

Prints the summary of a mnGLMM model.

## Usage

``` r
# S3 method for class 'mnGLMM'
summary(object, ...)
```

## Arguments

- object:

  An object of class "mnGLMM"

- ...:

  Additional arguments (currently unused).

## Value

`object`, invisibly. Called for the side effect of printing the
coefficient estimates from a fitted mnGLMM.

## Examples

``` r
data(story_pollen_matrix)
sample_idx <- which(rowSums(story_pollen_matrix) != 0)
n <- ncol(story_pollen_matrix)
glmm_mod <- mnGLMM(Y = story_pollen_matrix[sample_idx, ], V.fixed = diag(n))
summary(glmm_mod)
#> 
#> Call: mnGLMM with Tmax =  93  n =  5 
#> 
#> logLik = 449.3075,  AIC = -888.6149 [df = 5]
#> 
#> dispersion parameter =  1
#> 
#> Fitted Coefficients with approximate se
#>                                     Coef.         se         t            P
#> (intercept).hardwood           0.09345376 0.07082625  1.319479 1.870090e-01
#> (intercept).Fagus grandifolia -0.31836613 0.07195283 -4.424650 9.659855e-06
#> (intercept).Ulmus             -0.39757671 0.07161670 -5.551453 2.833052e-08
#> (intercept).Quercus            0.84814218 0.07232763 11.726393 9.335258e-32
#> 
#> 
#> Overall model
#> 
#> B =  
#>      [,1]       [,2]       [,3]       [,4]      [,5]
#> [1,]    0 0.09345376 -0.3183661 -0.3975767 0.8481422
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
