# Story Charcoal Matrix

A matrix of charcoal accumulation rates from Story Lake, Indiana. This
matrix has been binned at a 50-year resolution. Values are on the
original measurement scale; scale them (e.g. with
[`scale()`](https://rdrr.io/r/base/scale.html)) before passing them to
[`mnGLMM()`](https://quinnasena.github.io/multinomialTS/reference/mnGLMM.md)
or
[`mnTS()`](https://quinnasena.github.io/multinomialTS/reference/mnTS.md).

## Usage

``` r
story_char_matrix
```

## Format

A `1` column matrix with `160` rows (observations).

## Source

For details see vignette

## Examples

``` r
data(story_char_matrix)
head(story_char_matrix)
#>       char_acc
#> [1,]  3.428030
#> [2,]  4.002948
#> [3,] 11.145338
#> [4,]  4.060201
#> [5,]  5.292200
#> [6,]  8.241451
```
