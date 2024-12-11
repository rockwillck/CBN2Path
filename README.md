
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CBNpath

<!-- badges: start -->
<!-- badges: end -->

The goal of CBNpath is to provide an R interface to the CBN family of
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

You can install the development version of `CBNpath` like so:

``` r
devtools::install_github("rockwillck/rCBN")
```

## Example

`CBNpath` provides an interface to CT-CBN, H-CBN, and B-CBN:

``` r
library(CBNpath)
bc = Spock$new(
  poset = read_poset("inst/extdata/BC")$sets,
  numMutations = read_poset("inst/extdata/BC")$mutations,
  patternOrLambda = read_pattern("inst/extdata/BC")
)
```

``` r
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
#> 
#> $row
#>         Poset           Eps         Alpha       Loglike      lambda_s 
#>     0.0000000     1.0000000     1.0000000 -4731.7300000     1.0000000 
#>      lambda_1      lambda_2      lambda_3      lambda_4      lambda_5 
#>     1.2637990     0.2796195     0.3760628     0.8210815     0.3193844 
#>      lambda_6      lambda_7      lambda_8      lambda_9     lambda_10 
#>     0.3212656     0.3178370     0.3777232     0.5241736     0.4419216
```

``` r
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

``` r
bcbn()
#> [1] "chain: 1"
#> [1] 0
#> [1] "chain: 2"
#> [1] 0
#> [1] "chain: 3"
#> [1] 0
#> [1] "chain: 4"
#> [1] 0
#>        V1                V2               V3                  V4         
#>  Min.   :0.02107   Min.   :0.1352   Min.   :0.0003081   Min.   :0.02569  
#>  1st Qu.:0.66748   1st Qu.:0.7637   1st Qu.:0.4366964   1st Qu.:0.10239  
#>  Median :0.80378   Median :0.8723   Median :0.5969084   Median :0.12979  
#>  Mean   :0.77148   Mean   :0.8389   Mean   :0.5871786   Mean   :0.13266  
#>  3rd Qu.:0.90473   3rd Qu.:0.9476   3rd Qu.:0.7499415   3rd Qu.:0.15900  
#>  Max.   :0.99999   Max.   :1.0000   Max.   :0.9999900   Max.   :0.31579  
#>        V5         
#>  Min.   :-12.447  
#>  1st Qu.: -8.181  
#>  Median : -7.755  
#>  Mean   : -7.882  
#>  3rd Qu.: -7.453  
#>  Max.   : -6.991  
#>        V1                V2                V3                 V4         
#>  Min.   :0.06648   Min.   :0.04497   Min.   :0.003555   Min.   :0.02746  
#>  1st Qu.:0.66208   1st Qu.:0.76233   1st Qu.:0.439097   1st Qu.:0.10125  
#>  Median :0.79747   Median :0.87339   Median :0.599634   Median :0.12707  
#>  Mean   :0.76753   Mean   :0.83623   Mean   :0.589607   Mean   :0.13150  
#>  3rd Qu.:0.90465   3rd Qu.:0.94537   3rd Qu.:0.749613   3rd Qu.:0.15668  
#>  Max.   :0.99999   Max.   :0.99996   Max.   :0.999870   Max.   :0.32878  
#>        V5         
#>  Min.   :-17.714  
#>  1st Qu.: -8.185  
#>  Median : -7.759  
#>  Mean   : -7.887  
#>  3rd Qu.: -7.453  
#>  Max.   : -6.979  
#>        V1                V2                V3                 V4         
#>  Min.   :0.02013   Min.   :0.09462   Min.   :0.001363   Min.   :0.02796  
#>  1st Qu.:0.66826   1st Qu.:0.75962   1st Qu.:0.431286   1st Qu.:0.10154  
#>  Median :0.80737   Median :0.87060   Median :0.588016   Median :0.12879  
#>  Mean   :0.77365   Mean   :0.83551   Mean   :0.581040   Mean   :0.13140  
#>  3rd Qu.:0.91063   3rd Qu.:0.94250   3rd Qu.:0.740201   3rd Qu.:0.15698  
#>  Max.   :0.99997   Max.   :0.99997   Max.   :0.999796   Max.   :0.35379  
#>        V5         
#>  Min.   :-12.807  
#>  1st Qu.: -8.173  
#>  Median : -7.758  
#>  Mean   : -7.883  
#>  3rd Qu.: -7.458  
#>  Max.   : -6.972  
#>        V1                V2               V3                  V4         
#>  Min.   :0.05784   Min.   :0.1464   Min.   :0.0006188   Min.   :0.02212  
#>  1st Qu.:0.65974   1st Qu.:0.7686   1st Qu.:0.4360143   1st Qu.:0.10171  
#>  Median :0.79969   Median :0.8752   Median :0.5911433   Median :0.12858  
#>  Mean   :0.76814   Mean   :0.8418   Mean   :0.5852267   Mean   :0.13204  
#>  3rd Qu.:0.90621   3rd Qu.:0.9473   3rd Qu.:0.7445760   3rd Qu.:0.15728  
#>  Max.   :0.99996   Max.   :1.0000   Max.   :0.9998063   Max.   :0.33629  
#>        V5         
#>  Min.   :-12.383  
#>  1st Qu.: -8.173  
#>  Median : -7.761  
#>  Mean   : -7.880  
#>  3rd Qu.: -7.453  
#>  Max.   : -6.980  
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
#> [1] "finished run: 1"
#> [1] "chain: 1"
#> [1] 0
#> [1] "chain: 2"
#> [1] 0
#> [1] "chain: 3"
#> [1] 0
#> [1] "chain: 4"
#> [1] 0
#>        V1                  V2                V3                 V4         
#>  Min.   :0.0004549   Min.   :0.02496   Min.   :0.003726   Min.   :0.02906  
#>  1st Qu.:0.6678211   1st Qu.:0.76538   1st Qu.:0.439366   1st Qu.:0.10397  
#>  Median :0.7968230   Median :0.87450   Median :0.595538   Median :0.12979  
#>  Mean   :0.7681791   Mean   :0.83799   Mean   :0.587078   Mean   :0.13264  
#>  3rd Qu.:0.8997379   3rd Qu.:0.94537   3rd Qu.:0.744427   3rd Qu.:0.15874  
#>  Max.   :0.9999542   Max.   :0.99999   Max.   :0.999981   Max.   :0.31917  
#>        V5         
#>  Min.   :-12.907  
#>  1st Qu.: -8.147  
#>  Median : -7.734  
#>  Mean   : -7.861  
#>  3rd Qu.: -7.445  
#>  Max.   : -6.977  
#>        V1                V2                V3                 V4         
#>  Min.   :0.02725   Min.   :0.09019   Min.   :0.002087   Min.   :0.02721  
#>  1st Qu.:0.67057   1st Qu.:0.76377   1st Qu.:0.433980   1st Qu.:0.10216  
#>  Median :0.80489   Median :0.87300   Median :0.589834   Median :0.12737  
#>  Mean   :0.77276   Mean   :0.83794   Mean   :0.585187   Mean   :0.13264  
#>  3rd Qu.:0.90768   3rd Qu.:0.94623   3rd Qu.:0.744392   3rd Qu.:0.15690  
#>  Max.   :0.99996   Max.   :0.99998   Max.   :0.999769   Max.   :0.33502  
#>        V5         
#>  Min.   :-12.395  
#>  1st Qu.: -8.171  
#>  Median : -7.752  
#>  Mean   : -7.879  
#>  3rd Qu.: -7.454  
#>  Max.   : -6.974  
#>        V1                 V2                V3                 V4         
#>  Min.   :0.001597   Min.   :0.01512   Min.   :0.001117   Min.   :0.02424  
#>  1st Qu.:0.664525   1st Qu.:0.76598   1st Qu.:0.439307   1st Qu.:0.10092  
#>  Median :0.796754   Median :0.87326   Median :0.589886   Median :0.12759  
#>  Mean   :0.770036   Mean   :0.83755   Mean   :0.586562   Mean   :0.13138  
#>  3rd Qu.:0.906064   3rd Qu.:0.94577   3rd Qu.:0.743802   3rd Qu.:0.15740  
#>  Max.   :0.999990   Max.   :0.99998   Max.   :0.999969   Max.   :0.33868  
#>        V5         
#>  Min.   :-11.902  
#>  1st Qu.: -8.164  
#>  Median : -7.753  
#>  Mean   : -7.876  
#>  3rd Qu.: -7.453  
#>  Max.   : -6.978  
#>        V1                V2               V3                  V4         
#>  Min.   :0.05311   Min.   :0.1489   Min.   :0.0002457   Min.   :0.01653  
#>  1st Qu.:0.67016   1st Qu.:0.7648   1st Qu.:0.4336892   1st Qu.:0.10370  
#>  Median :0.80202   Median :0.8742   Median :0.5918887   Median :0.12919  
#>  Mean   :0.77249   Mean   :0.8397   Mean   :0.5842680   Mean   :0.13283  
#>  3rd Qu.:0.90597   3rd Qu.:0.9470   3rd Qu.:0.7456341   3rd Qu.:0.15897  
#>  Max.   :1.00000   Max.   :1.0000   Max.   :0.9998586   Max.   :0.35206  
#>        V5         
#>  Min.   :-12.984  
#>  1st Qu.: -8.157  
#>  Median : -7.749  
#>  Mean   : -7.868  
#>  3rd Qu.: -7.452  
#>  Max.   : -6.980  
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
#>         [,1]    [,2]    [,3]
#> [1,] 0.00000 0.12418 0.33958
#> [2,] 0.26901 0.00000 0.42452
#> [3,] 0.04919 0.02925 0.00000
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
