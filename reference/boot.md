# Bootstrap a model (deprecated)

Deprecated in favour of
[`bootstrap()`](https://quinnasena.github.io/multinomialTS/reference/bootstrap.md),
which does the same thing under a name that does not collide with
[`boot::boot()`](https://rdrr.io/pkg/boot/man/boot.html) from the
recommended boot package. This wrapper warns and forwards, and will be
removed in a future release.

## Usage

``` r
boot(object, ...)
```

## Arguments

- object:

  An object of class `mnTS`.

- ...:

  Passed to
  [`bootstrap()`](https://quinnasena.github.io/multinomialTS/reference/bootstrap.md).

## Value

The value returned by
[`bootstrap()`](https://quinnasena.github.io/multinomialTS/reference/bootstrap.md).

## See also

[`bootstrap()`](https://quinnasena.github.io/multinomialTS/reference/bootstrap.md),
which replaces this function.
