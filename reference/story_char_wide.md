# Story Charcoal Wide

A ~10,000-year data set from Story Lake, Indiana, of modelled ages,
sample depths, and charcoal accumulation rate from a sediment core
sample.

## Usage

``` r
story_char_wide
```

## Format

A tibble with `736` rows and `3` columns.

## Source

For details see https://doi.org/10.1111/1365-2745.14289

## Examples

``` r
data(story_char_wide)
str(story_char_wide)
#> tibble [736 × 3] (S3: tbl_df/tbl/data.frame)
#>  $ age     : num [1:736] -66.2 -63.4 -61.5 -59.8 -58.5 ...
#>  $ depth   : num [1:736] 0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 ...
#>  $ char_acc: num [1:736] 0 15.8 40.3 21.8 37.3 ...
```
