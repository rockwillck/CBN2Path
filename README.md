
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CBN2Path

<!-- badges: start -->
<!-- badges: end -->

The goal of CBN2Path is to provide an R interface to the CBN family of
functions developed at ETH-Zurich.

## Dependencies

`gsl` can be installed using:

``` bash
brew install gsl
```

If `gsl` was installed using any method other Homebrew, uninstall `gsl`
before reinstalling using Homebrew.

If Homebrew is not installed, follow the instructions at
<https://brew.sh> to install.

Make sure to restart R before proceeding.

## Installation

You can install the development version of `CBN2Path` like so:

``` r
devtools::install_github("rockwillck/CBN2Path")
```

## Example

`CBN2Path` provides an interface to CT-CBN, H-CBN, and B-CBN:

``` r
library(CBN2Path)
bc <- Spock$new(
    poset = read_poset("inst/extdata/BC")$sets,
    numMutations = read_poset("inst/extdata/BC")$mutations,
    genotypeMatrix = read_pattern("inst/extdata/BC")
)
```

``` r
ctcbn(bc)
#> [[1]]
#> [[1]]$lambda
#>           [,1]
#>  [1,] 1.000000
#>  [2,] 1.338916
#>  [3,] 0.642515
#>  [4,] 1.093539
#>  [5,] 0.973560
#>  [6,] 1.424189
#>  [7,] 0.981581
#>  [8,] 0.539870
#>  [9,] 0.438919
#> [10,] 0.787104
#> [11,] 1.452059
#> 
#> [[1]]$summary
#>         Poset           Eps         Alpha       Loglike      lambda_s 
#>     0.0000000     0.0000000     0.4871480 -5208.5600000     1.0000000 
#>      lambda_1      lambda_2      lambda_3      lambda_4      lambda_5 
#>     1.3389157     0.6425153     1.0935394     0.9735600     1.4241889 
#>      lambda_6      lambda_7      lambda_8      lambda_9     lambda_10 
#>     0.9815809     0.5398705     0.4389191     0.7871042     1.4520586
```

``` r
hcbn(bc)
#> [[1]]
#> [[1]]$poset
#> [[1]]$poset$mutations
#> [1] 10
#> 
#> [[1]]$poset$sets
#>       [,1] [,2]
#>  [1,]    1    2
#>  [2,]    1    3
#>  [3,]    1    8
#>  [4,]    3    5
#>  [5,]    4    2
#>  [6,]    4    3
#>  [7,]    9    7
#>  [8,]    9   10
#>  [9,]   10    6
#> 
#> 
#> [[1]]$lambda
#>           [,1]
#>  [1,] 1.000000
#>  [2,] 1.658991
#>  [3,] 0.473723
#>  [4,] 1.219269
#>  [5,] 0.889385
#>  [6,] 2.584040
#>  [7,] 1.718128
#>  [8,] 0.750927
#>  [9,] 0.472000
#> [10,] 0.463375
#> [11,] 2.967137
#> 
#> [[1]]$summary
#>         Poset           Eps         Alpha       Loglike      lambda_s 
#>     0.0000000     0.1505320     0.4871480 -4677.6000000     1.0000000 
#>      lambda_1      lambda_2      lambda_3      lambda_4      lambda_5 
#>     1.6589908     0.4737230     1.2192686     0.8893852     2.5840398 
#>      lambda_6      lambda_7      lambda_8      lambda_9     lambda_10 
#>     1.7181281     0.7509271     0.4720001     0.4633753     2.9671369
```

``` r
bcbn()
#> Loading required package: rBCBN
#> [1] "chain: 1"
#> [1] 0
#> [1] "chain: 2"
#> [1] 0
#> [1] "chain: 3"
#> [1] 0
#> [1] "chain: 4"
#> [1] 0
#>        V1               V2                 V3                  V4         
#>  Min.   :0.1186   Min.   :0.001451   Min.   :0.0000422   Min.   :0.02232  
#>  1st Qu.:0.8107   1st Qu.:0.474243   1st Qu.:0.2653763   1st Qu.:0.09170  
#>  Median :0.8982   Median :0.605459   Median :0.4285352   Median :0.11730  
#>  Mean   :0.8680   Mean   :0.602572   Mean   :0.4447468   Mean   :0.12095  
#>  3rd Qu.:0.9570   3rd Qu.:0.740179   3rd Qu.:0.6104647   3rd Qu.:0.14661  
#>  Max.   :1.0000   Max.   :0.999531   Max.   :0.9998044   Max.   :0.32046  
#>        V5         
#>  Min.   :-11.612  
#>  1st Qu.: -7.817  
#>  Median : -7.393  
#>  Mean   : -7.498  
#>  3rd Qu.: -7.055  
#>  Max.   : -6.533  
#>        V1                V2                V3                  V4         
#>  Min.   :0.08362   Min.   :0.01064   Min.   :0.0002166   Min.   :0.01639  
#>  1st Qu.:0.80999   1st Qu.:0.47424   1st Qu.:0.2648785   1st Qu.:0.09092  
#>  Median :0.89912   Median :0.60581   Median :0.4296266   Median :0.11753  
#>  Mean   :0.86778   Mean   :0.60295   Mean   :0.4454813   Mean   :0.12157  
#>  3rd Qu.:0.95766   3rd Qu.:0.73784   3rd Qu.:0.6137740   3rd Qu.:0.14745  
#>  Max.   :0.99999   Max.   :0.99984   Max.   :0.9995950   Max.   :0.33986  
#>        V5         
#>  Min.   :-11.643  
#>  1st Qu.: -7.822  
#>  Median : -7.395  
#>  Mean   : -7.504  
#>  3rd Qu.: -7.057  
#>  Max.   : -6.537  
#>        V1               V2                 V3                  V4         
#>  Min.   :0.2056   Min.   :0.005852   Min.   :0.0000238   Min.   :0.02032  
#>  1st Qu.:0.8068   1st Qu.:0.467740   1st Qu.:0.2715890   1st Qu.:0.09241  
#>  Median :0.8982   Median :0.606022   Median :0.4391039   Median :0.11743  
#>  Mean   :0.8655   Mean   :0.602039   Mean   :0.4510079   Mean   :0.12114  
#>  3rd Qu.:0.9539   3rd Qu.:0.738240   3rd Qu.:0.6198538   3rd Qu.:0.14639  
#>  Max.   :1.0000   Max.   :0.999903   Max.   :0.9997733   Max.   :0.34131  
#>        V5         
#>  Min.   :-11.460  
#>  1st Qu.: -7.839  
#>  Median : -7.405  
#>  Mean   : -7.517  
#>  3rd Qu.: -7.078  
#>  Max.   : -6.536  
#>        V1               V2                  V3                  V4         
#>  Min.   :0.1201   Min.   :0.0005872   Min.   :0.0000019   Min.   :0.01837  
#>  1st Qu.:0.8088   1st Qu.:0.4729241   1st Qu.:0.2685715   1st Qu.:0.09015  
#>  Median :0.8983   Median :0.6074857   Median :0.4330691   Median :0.11704  
#>  Mean   :0.8672   Mean   :0.6025209   Mean   :0.4475112   Mean   :0.12029  
#>  3rd Qu.:0.9556   3rd Qu.:0.7373943   3rd Qu.:0.6143208   3rd Qu.:0.14550  
#>  Max.   :1.0000   Max.   :0.9997110   Max.   :0.9999610   Max.   :0.31443  
#>        V5         
#>  Min.   :-11.978  
#>  1st Qu.: -7.819  
#>  Median : -7.388  
#>  Mean   : -7.498  
#>  3rd Qu.: -7.054  
#>  Max.   : -6.537  
#> [1] "Criterion: 1.00054325997706"
#> Potential scale reduction factors:
#> 
#>      Point est. Upper C.I.
#> [1,]          1          1
#> [2,]          1          1
#> [3,]          1          1
#> [4,]          1          1
#> [5,]          1          1
#> 
#> Multivariate psrf
#> 
#> 1
#> [1] "##########################################"
#>         [,1]    [,2]    [,3]
#> [1,] 0.00000 0.22664 0.32146
#> [2,] 0.00947 0.00000 0.70233
#> [3,] 0.00012 0.05358 0.00000
#>         [,1]    [,2]    [,3]
#> [1,] 0.00000 0.22664 0.32146
#> [2,] 0.00947 0.00000 0.70233
#> [3,] 0.00012 0.05358 0.00000
```

## CT-CBN

This project is based on and draws on code from CT-CBN from ETH
Zurich.  
CT-CBN. (n.d.). Department of Biosystems Science and Engineering –
Department of Biosystems Science and Engineering \| ETH Zurich.
<https://bsse.ethz.ch/cbg/software/ct-cbn.html>

## B-CBN

This project is based on and draws on code from B-CBN from ETH Zurich.  
bcbn. (n.d.). Department of Biosystems Science and Engineering –
Department of Biosystems Science and Engineering \| ETH Zurich.
<https://bsse.ethz.ch/cbg/software/bcbn.html>
