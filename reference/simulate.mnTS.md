# Simulate Data from a Multinomial State-Space Model

Generates simulated datasets from a fitted mnTS object.

## Usage

``` r
# S3 method for class 'mnTS'
simulate(object, ...)
```

## Arguments

- object:

  An object of class `"mnTS"`, as returned by the
  [`mnTS`](https://quinnasena.github.io/multinomialTS/reference/mnTS.md)
  function.

- ...:

  Additional arguments passed to the simulation routine.

## Value

A simulated dataset in the same structure as the response variable `Y`
used to fit the model.
