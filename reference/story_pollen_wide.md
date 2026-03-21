# Story Pollen Wide

A ~10,000-year data set from Story Lake, Indiana, of modelled ages,
sample depths, pollen taxa. Data are organised into a reference group of
"other" taxa, a functional group of hardwood taxa, and target taxa of
Fagus grandifolia, Ulmus spp., and Quercus spp.

## Usage

``` r
story_pollen_wide
```

## Format

A tibble with \`95\` rows (observations) and \`7\` columns.

## Source

For details see https://doi.org/10.1111/1365-2745.14289

## Examples

``` r
data(story_pollen_wide)
names(story_pollen_wide)
#> [1] "age"               "depth"             "other"            
#> [4] "hardwood"          "Fagus grandifolia" "Quercus"          
#> [7] "Ulmus"            
```
