# Story Pollen Matrix

A matrix of pollen taxa counts from Story Lake, Indiana. Data are
organised into a reference group of "other" taxa, a functional group of
hardwood taxa, and target taxa of Fagus grandifolia, Ulmus spp., and
Quercus spp. Data have been binned at a 50-year resolution with empty
bins filled with 0-values (i.e., rows of 0 are empty bins, see vignette
for details).

## Usage

``` r
story_pollen_matrix
```

## Format

A matrix with \`160\` rows (samples) and \`5\` columns (taxa).

## Source

For details see vignette

## Examples

``` r
data(story_pollen_matrix)
summary(story_pollen_matrix)
#>      other           hardwood      Fagus grandifolia     Ulmus      
#>  Min.   :  0.00   Min.   :  0.00   Min.   : 0.00     Min.   : 0.00  
#>  1st Qu.:  0.00   1st Qu.:  0.00   1st Qu.: 0.00     1st Qu.: 0.00  
#>  Median : 30.00   Median : 48.00   Median : 7.50     Median :16.50  
#>  Mean   : 26.94   Mean   : 35.66   Mean   :19.71     Mean   :18.87  
#>  3rd Qu.: 40.00   3rd Qu.: 62.25   3rd Qu.:39.00     3rd Qu.:29.00  
#>  Max.   :291.00   Max.   :161.00   Max.   :98.00     Max.   :84.00  
#>     Quercus      
#>  Min.   :  0.00  
#>  1st Qu.:  0.00  
#>  Median : 91.50  
#>  Mean   : 70.31  
#>  3rd Qu.:121.25  
#>  Max.   :248.00  
```
