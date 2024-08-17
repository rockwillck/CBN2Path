
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

rCBN provides an interface to CT-CBN, H-CBN, and B-CBN:

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
#>        V1                V2                 V3                  V4         
#>  Min.   :0.06924   Min.   :0.009552   Min.   :0.0000469   Min.   :0.02316  
#>  1st Qu.:0.72912   1st Qu.:0.561202   1st Qu.:0.3346137   1st Qu.:0.11230  
#>  Median :0.84907   Median :0.715177   Median :0.4959095   Median :0.13940  
#>  Mean   :0.81317   Mean   :0.691535   Mean   :0.4995369   Mean   :0.14373  
#>  3rd Qu.:0.93380   3rd Qu.:0.847162   3rd Qu.:0.6631499   3rd Qu.:0.17318  
#>  Max.   :0.99999   Max.   :0.999958   Max.   :0.9996534   Max.   :0.33584  
#>        V5         
#>  Min.   :-12.907  
#>  1st Qu.: -9.072  
#>  Median : -8.670  
#>  Mean   : -8.796  
#>  3rd Qu.: -8.393  
#>  Max.   : -7.878  
#>        V1               V2                 V3                  V4         
#>  Min.   :0.1173   Min.   :0.007604   Min.   :0.0000234   Min.   :0.02374  
#>  1st Qu.:0.7341   1st Qu.:0.554134   1st Qu.:0.3301480   1st Qu.:0.10788  
#>  Median :0.8565   Median :0.708195   Median :0.4923060   Median :0.13825  
#>  Mean   :0.8196   Mean   :0.686655   Mean   :0.4974030   Mean   :0.14238  
#>  3rd Qu.:0.9399   3rd Qu.:0.842289   3rd Qu.:0.6605243   3rd Qu.:0.17214  
#>  Max.   :1.0000   Max.   :0.999985   Max.   :0.9999242   Max.   :0.35538  
#>        V5         
#>  Min.   :-13.502  
#>  1st Qu.: -9.081  
#>  Median : -8.689  
#>  Mean   : -8.805  
#>  3rd Qu.: -8.401  
#>  Max.   : -7.872  
#>        V1                V2               V3                  V4         
#>  Min.   :0.03585   Min.   :0.0129   Min.   :0.0000611   Min.   :0.02855  
#>  1st Qu.:0.72844   1st Qu.:0.5541   1st Qu.:0.3353986   1st Qu.:0.11322  
#>  Median :0.85231   Median :0.7105   Median :0.4984755   Median :0.14253  
#>  Mean   :0.81508   Mean   :0.6886   Mean   :0.5003764   Mean   :0.14493  
#>  3rd Qu.:0.93815   3rd Qu.:0.8474   3rd Qu.:0.6617609   3rd Qu.:0.17367  
#>  Max.   :0.99993   Max.   :1.0000   Max.   :0.9999534   Max.   :0.33031  
#>        V5         
#>  Min.   :-13.124  
#>  1st Qu.: -9.054  
#>  Median : -8.660  
#>  Mean   : -8.786  
#>  3rd Qu.: -8.384  
#>  Max.   : -7.872  
#>        V1                V2                V3                  V4         
#>  Min.   :0.03513   Min.   :0.01018   Min.   :0.0000762   Min.   :0.01977  
#>  1st Qu.:0.73121   1st Qu.:0.55838   1st Qu.:0.3324746   1st Qu.:0.11106  
#>  Median :0.84939   Median :0.71454   Median :0.4897298   Median :0.13883  
#>  Mean   :0.81454   Mean   :0.69074   Mean   :0.4960991   Mean   :0.14295  
#>  3rd Qu.:0.93572   3rd Qu.:0.84598   3rd Qu.:0.6575516   3rd Qu.:0.17244  
#>  Max.   :0.99992   Max.   :1.00000   Max.   :0.9997885   Max.   :0.34139  
#>        V5         
#>  Min.   :-13.425  
#>  1st Qu.: -9.058  
#>  Median : -8.659  
#>  Mean   : -8.786  
#>  3rd Qu.: -8.385  
#>  Max.   : -7.871  
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
#>        V1                V2                 V3                  V4         
#>  Min.   :0.01528   Min.   :0.005781   Min.   :0.0000415   Min.   :0.02489  
#>  1st Qu.:0.73274   1st Qu.:0.557696   1st Qu.:0.3340968   1st Qu.:0.11092  
#>  Median :0.85248   Median :0.713161   Median :0.4956646   Median :0.14059  
#>  Mean   :0.81438   Mean   :0.690483   Mean   :0.5004620   Mean   :0.14377  
#>  3rd Qu.:0.93410   3rd Qu.:0.848504   3rd Qu.:0.6664166   3rd Qu.:0.17259  
#>  Max.   :0.99991   Max.   :0.999871   Max.   :0.9999581   Max.   :0.35988  
#>        V5         
#>  Min.   :-12.625  
#>  1st Qu.: -9.063  
#>  Median : -8.673  
#>  Mean   : -8.792  
#>  3rd Qu.: -8.396  
#>  Max.   : -7.868  
#>        V1               V2                 V3                  V4         
#>  Min.   :0.0642   Min.   :0.001336   Min.   :0.0000252   Min.   :0.02911  
#>  1st Qu.:0.7312   1st Qu.:0.561069   1st Qu.:0.3301559   1st Qu.:0.11186  
#>  Median :0.8524   Median :0.711602   Median :0.4927934   Median :0.14005  
#>  Mean   :0.8144   Mean   :0.691280   Mean   :0.4989716   Mean   :0.14358  
#>  3rd Qu.:0.9337   3rd Qu.:0.850244   3rd Qu.:0.6652211   3rd Qu.:0.17107  
#>  Max.   :1.0000   Max.   :0.999970   Max.   :0.9999392   Max.   :0.33184  
#>        V5         
#>  Min.   :-13.456  
#>  1st Qu.: -9.046  
#>  Median : -8.660  
#>  Mean   : -8.781  
#>  3rd Qu.: -8.389  
#>  Max.   : -7.865  
#>        V1                 V2                V3                  V4         
#>  Min.   :0.008725   Min.   :0.01276   Min.   :0.0003281   Min.   :0.02359  
#>  1st Qu.:0.733368   1st Qu.:0.55866   1st Qu.:0.3373365   1st Qu.:0.10949  
#>  Median :0.853101   Median :0.71389   Median :0.4959790   Median :0.13803  
#>  Mean   :0.815844   Mean   :0.69139   Mean   :0.5006825   Mean   :0.14257  
#>  3rd Qu.:0.935814   3rd Qu.:0.84666   3rd Qu.:0.6630031   3rd Qu.:0.17037  
#>  Max.   :0.999997   Max.   :0.99994   Max.   :0.9999512   Max.   :0.37880  
#>        V5         
#>  Min.   :-13.979  
#>  1st Qu.: -9.071  
#>  Median : -8.673  
#>  Mean   : -8.792  
#>  3rd Qu.: -8.394  
#>  Max.   : -7.877  
#>        V1                 V2                  V3                  V4         
#>  Min.   :0.002852   Min.   :0.0001984   Min.   :0.0001913   Min.   :0.03181  
#>  1st Qu.:0.732462   1st Qu.:0.5570014   1st Qu.:0.3347501   1st Qu.:0.11229  
#>  Median :0.853659   Median :0.7118792   Median :0.4959840   Median :0.14029  
#>  Mean   :0.818568   Mean   :0.6901617   Mean   :0.4988545   Mean   :0.14443  
#>  3rd Qu.:0.939435   3rd Qu.:0.8467971   3rd Qu.:0.6623643   3rd Qu.:0.17301  
#>  Max.   :0.999969   Max.   :0.9999319   Max.   :0.9999685   Max.   :0.35606  
#>        V5         
#>  Min.   :-12.924  
#>  1st Qu.: -9.060  
#>  Median : -8.659  
#>  Mean   : -8.783  
#>  3rd Qu.: -8.387  
#>  Max.   : -7.874  
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
#> [1,] 0.00000 0.34100 0.28146
#> [2,] 0.10297 0.00000 0.34569
#> [3,] 0.01160 0.08407 0.00000
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
