# Bootstrap method for mnTS objects

Provides bootstrapped confidence intervals for an object of class
`mnTS`.

## Usage

``` r
# S3 method for class 'mnTS'
boot(mod, reps, dispersion.fixed = 1, maxit.optim = 1e+05, ...)
```

## Arguments

- reps:

  Integer. Number of bootstrap replications.

- dispersion.fixed:

  Numeric. Dispersion parameter (default is 1).

- maxit.optim:

  Integer. Max iterations for optimization (default 100000).

- ...:

  Additional arguments (not currently used).

- object:

  An object of class `mnTS`.
