
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rCBN

<!-- badges: start -->
<!-- badges: end -->

The goal of rCBN is to provide an R interface to the CBN family of
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

You can install the development version of `rCBN` like so:

``` r
devtools::install_github("rockwillck/rCBN")
```

## Example

`rCBN` provides an interface to CT-CBN, H-CBN, and B-CBN:

``` r
library(rCBN)
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
#>        V1               V2                 V3                  V4         
#>  Min.   :0.1922   Min.   :0.009176   Min.   :0.0008584   Min.   :0.01480  
#>  1st Qu.:0.8572   1st Qu.:0.479691   1st Qu.:0.4302900   1st Qu.:0.07855  
#>  Median :0.9243   Median :0.611178   Median :0.5993672   Median :0.10232  
#>  Mean   :0.8987   Mean   :0.613320   Mean   :0.5912503   Mean   :0.10598  
#>  3rd Qu.:0.9679   3rd Qu.:0.749452   3rd Qu.:0.7617102   3rd Qu.:0.12832  
#>  Max.   :1.0000   Max.   :0.999789   Max.   :0.9999754   Max.   :0.28584  
#>        V5         
#>  Min.   :-13.347  
#>  1st Qu.: -6.982  
#>  Median : -6.495  
#>  Mean   : -6.608  
#>  3rd Qu.: -6.104  
#>  Max.   : -5.543  
#>        V1               V2                V3                  V4         
#>  Min.   :0.1512   Min.   :0.01631   Min.   :0.0005482   Min.   :0.01643  
#>  1st Qu.:0.8557   1st Qu.:0.48167   1st Qu.:0.4333057   1st Qu.:0.07895  
#>  Median :0.9263   Median :0.61505   Median :0.6049773   Median :0.10207  
#>  Mean   :0.8968   Mean   :0.61675   Mean   :0.5938923   Mean   :0.10626  
#>  3rd Qu.:0.9684   3rd Qu.:0.75131   3rd Qu.:0.7651331   3rd Qu.:0.12951  
#>  Max.   :1.0000   Max.   :0.99992   Max.   :0.9999885   Max.   :0.30495  
#>        V5         
#>  Min.   :-11.379  
#>  1st Qu.: -6.990  
#>  Median : -6.480  
#>  Mean   : -6.604  
#>  3rd Qu.: -6.097  
#>  Max.   : -5.543  
#>        V1               V2                V3                  V4         
#>  Min.   :0.1696   Min.   :0.01848   Min.   :0.0007269   Min.   :0.01325  
#>  1st Qu.:0.8613   1st Qu.:0.47856   1st Qu.:0.4348357   1st Qu.:0.07901  
#>  Median :0.9299   Median :0.61133   Median :0.6002524   Median :0.10216  
#>  Mean   :0.9027   Mean   :0.61237   Mean   :0.5944681   Mean   :0.10572  
#>  3rd Qu.:0.9713   3rd Qu.:0.74588   3rd Qu.:0.7706866   3rd Qu.:0.12844  
#>  Max.   :1.0000   Max.   :0.99989   Max.   :0.9998931   Max.   :0.30845  
#>        V5         
#>  Min.   :-14.949  
#>  1st Qu.: -6.964  
#>  Median : -6.466  
#>  Mean   : -6.592  
#>  3rd Qu.: -6.080  
#>  Max.   : -5.549  
#>        V1               V2                V3                  V4         
#>  Min.   :0.2927   Min.   :0.01797   Min.   :0.0008176   Min.   :0.01526  
#>  1st Qu.:0.8580   1st Qu.:0.48188   1st Qu.:0.4362692   1st Qu.:0.07718  
#>  Median :0.9272   Median :0.61447   Median :0.6042587   Median :0.10037  
#>  Mean   :0.8999   Mean   :0.61346   Mean   :0.5958655   Mean   :0.10498  
#>  3rd Qu.:0.9684   3rd Qu.:0.74615   3rd Qu.:0.7692735   3rd Qu.:0.12766  
#>  Max.   :1.0000   Max.   :0.99968   Max.   :0.9998902   Max.   :0.33141  
#>        V5         
#>  Min.   :-11.536  
#>  1st Qu.: -6.960  
#>  Median : -6.467  
#>  Mean   : -6.594  
#>  3rd Qu.: -6.090  
#>  Max.   : -5.538  
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
#>        V1               V2               V3                  V4         
#>  Min.   :0.2909   Min.   :0.0013   Min.   :0.0002673   Min.   :0.01635  
#>  1st Qu.:0.8632   1st Qu.:0.4809   1st Qu.:0.4354553   1st Qu.:0.07783  
#>  Median :0.9281   Median :0.6106   Median :0.6006505   Median :0.10152  
#>  Mean   :0.9021   Mean   :0.6124   Mean   :0.5939053   Mean   :0.10453  
#>  3rd Qu.:0.9695   3rd Qu.:0.7453   3rd Qu.:0.7677360   3rd Qu.:0.12808  
#>  Max.   :1.0000   Max.   :0.9998   Max.   :0.9999861   Max.   :0.30394  
#>        V5         
#>  Min.   :-10.651  
#>  1st Qu.: -6.980  
#>  Median : -6.479  
#>  Mean   : -6.593  
#>  3rd Qu.: -6.086  
#>  Max.   : -5.544  
#>        V1               V2                 V3                  V4         
#>  Min.   :0.3143   Min.   :0.006877   Min.   :0.0000604   Min.   :0.01765  
#>  1st Qu.:0.8601   1st Qu.:0.481646   1st Qu.:0.4353213   1st Qu.:0.07895  
#>  Median :0.9289   Median :0.612846   Median :0.5991983   Median :0.10279  
#>  Mean   :0.9006   Mean   :0.612952   Mean   :0.5923672   Mean   :0.10618  
#>  3rd Qu.:0.9696   3rd Qu.:0.746155   3rd Qu.:0.7615328   3rd Qu.:0.12909  
#>  Max.   :1.0000   Max.   :0.999967   Max.   :0.9999687   Max.   :0.31746  
#>        V5         
#>  Min.   :-10.562  
#>  1st Qu.: -6.966  
#>  Median : -6.477  
#>  Mean   : -6.596  
#>  3rd Qu.: -6.086  
#>  Max.   : -5.539  
#>        V1               V2                V3                 V4         
#>  Min.   :0.2559   Min.   :0.05184   Min.   :0.001182   Min.   :0.01725  
#>  1st Qu.:0.8594   1st Qu.:0.48037   1st Qu.:0.424779   1st Qu.:0.07881  
#>  Median :0.9287   Median :0.60999   Median :0.601786   Median :0.10245  
#>  Mean   :0.9003   Mean   :0.61375   Mean   :0.593196   Mean   :0.10555  
#>  3rd Qu.:0.9706   3rd Qu.:0.74640   3rd Qu.:0.773276   3rd Qu.:0.12921  
#>  Max.   :1.0000   Max.   :0.99996   Max.   :0.999935   Max.   :0.28158  
#>        V5         
#>  Min.   :-11.233  
#>  1st Qu.: -6.962  
#>  Median : -6.464  
#>  Mean   : -6.586  
#>  3rd Qu.: -6.087  
#>  Max.   : -5.541  
#>        V1               V2                 V3                  V4         
#>  Min.   :0.3535   Min.   :0.007728   Min.   :0.0000775   Min.   :0.01288  
#>  1st Qu.:0.8608   1st Qu.:0.480650   1st Qu.:0.4300956   1st Qu.:0.08061  
#>  Median :0.9281   Median :0.611575   Median :0.6015193   Median :0.10328  
#>  Mean   :0.9016   Mean   :0.613924   Mean   :0.5931135   Mean   :0.10716  
#>  3rd Qu.:0.9711   3rd Qu.:0.748698   3rd Qu.:0.7697697   3rd Qu.:0.12982  
#>  Max.   :1.0000   Max.   :0.999973   Max.   :0.9999871   Max.   :0.31058  
#>        V5         
#>  Min.   :-11.224  
#>  1st Qu.: -6.977  
#>  Median : -6.474  
#>  Mean   : -6.596  
#>  3rd Qu.: -6.089  
#>  Max.   : -5.548  
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
#> [1,] 0.00000 0.29721 0.35529
#> [2,] 0.00319 0.00000 0.72623
#> [3,] 0.00018 0.13979 0.00000
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
