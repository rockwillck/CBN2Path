
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rCBN

<!-- badges: start -->
<!-- badges: end -->

The goal of rCBN is to provide an R interface to the CBN family of
functions developed at ETH-Zurich.

## Installation

You can install the development version of rCBN like so:

``` r
devtools::install_github("rockwillck/rCBN")
```

## Example

rCBN provides an interface to both CT-CBN and H-CBN:

``` r
library(rCBN)

bc = Spock$new(
  poset = read_poset("../PATH/TO/examples/BC")$sets,
  numMutations = read_poset("../PATH/TO/examples/BC")$mutations,
  patternOrLambda = read_pattern("../PATH/TO/examples/BC")
)
ctcbn(bc)
#> $lambda
#>           [,1]
#>  [1,] 1.000000
#>  [2,] 1.263799
#>  [3,] 0.279620
#>  [4,] 0.376063
#>  [5,] 0.821082
#>  [6,] 0.319384
#>  [7,] 0.321266
#>  [8,] 0.317837
#>  [9,] 0.377723
#> [10,] 0.524174
#> [11,] 0.441922

hcbn(bc)
#> $poset
#> $poset$mutations
#> [1] 10
#> 
#> $poset$sets
#>      [,1]
#> [1,]   NA
#> 
#> 
#> $lambda
#>           [,1]
#>  [1,] 1.000000
#>  [2,] 1.263857
#>  [3,] 0.279615
#>  [4,] 0.376053
#>  [5,] 0.821075
#>  [6,] 0.319394
#>  [7,] 0.321268
#>  [8,] 0.317814
#>  [9,] 0.377686
#> [10,] 0.524145
#> [11,] 0.441899
```
