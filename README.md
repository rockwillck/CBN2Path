
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CBN2Path

<!-- badges: start -->
<!-- badges: end -->

## Authors:
William Choi-Kim and Sayed-Rzgar Hosseini

## Abstract:



## Installation

Before installing the package make sure that you have installed `gsl` as follows:

``` bash
brew install gsl
```
Note that if `gsl` was installed using any method other Homebrew, uninstall `gsl`
before reinstalling using Homebrew.
Furthermore, note that if Homebrew is not installed, follow the instructions at <https://brew.sh> to install.
Finally, make sure to restart R before proceeding.


Then you can install the development version of `CBN2Path` like so:

``` r
devtools::install_github("rockwillck/CBN2Path")
```

## Example

`CBN2Path` provides an interface to CT-CBN, H-CBN, and B-CBN:

``` r
library(CBN2Path)
bc = Spock$new(
  poset = read_poset("inst/extdata/BC")$sets,
  numMutations = read_poset("inst/extdata/BC")$mutations,
  genotypeMatrix = read_pattern("inst/extdata/BC")
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
#>        V1                V2                V3                  V4         
#>  Min.   :0.03865   Min.   :0.02386   Min.   :0.0003991   Min.   :0.02091  
#>  1st Qu.:0.68377   1st Qu.:0.66608   1st Qu.:0.5619727   1st Qu.:0.10244  
#>  Median :0.80958   Median :0.79778   Median :0.7077167   Median :0.13052  
#>  Mean   :0.78119   Mean   :0.76865   Mean   :0.6881920   Mean   :0.13402  
#>  3rd Qu.:0.90786   3rd Qu.:0.89771   3rd Qu.:0.8333036   3rd Qu.:0.16123  
#>  Max.   :0.99990   Max.   :0.99994   Max.   :0.9999935   Max.   :0.36736  
#>        V5         
#>  Min.   :-14.772  
#>  1st Qu.: -8.514  
#>  Median : -8.119  
#>  Mean   : -8.242  
#>  3rd Qu.: -7.835  
#>  Max.   : -7.272  
#>        V1                 V2                V3               V4         
#>  Min.   :0.006523   Min.   :0.03657   Min.   :0.0120   Min.   :0.02162  
#>  1st Qu.:0.689446   1st Qu.:0.66716   1st Qu.:0.5595   1st Qu.:0.10152  
#>  Median :0.812706   Median :0.79332   Median :0.7083   Median :0.12922  
#>  Mean   :0.783377   Mean   :0.76642   Mean   :0.6888   Mean   :0.13265  
#>  3rd Qu.:0.911142   3rd Qu.:0.89666   3rd Qu.:0.8371   3rd Qu.:0.16038  
#>  Max.   :0.999869   Max.   :1.00000   Max.   :1.0000   Max.   :0.33225  
#>        V5         
#>  Min.   :-12.986  
#>  1st Qu.: -8.523  
#>  Median : -8.119  
#>  Mean   : -8.244  
#>  3rd Qu.: -7.836  
#>  Max.   : -7.289  
#>        V1                 V2                V3                V4         
#>  Min.   :0.002269   Min.   :0.05515   Min.   :0.01467   Min.   :0.01696  
#>  1st Qu.:0.683829   1st Qu.:0.67074   1st Qu.:0.56625   1st Qu.:0.10369  
#>  Median :0.807566   Median :0.79798   Median :0.70682   Median :0.12960  
#>  Mean   :0.778592   Mean   :0.77024   Mean   :0.68865   Mean   :0.13428  
#>  3rd Qu.:0.906157   3rd Qu.:0.89719   3rd Qu.:0.83484   3rd Qu.:0.16105  
#>  Max.   :0.999963   Max.   :0.99998   Max.   :0.99982   Max.   :0.34033  
#>        V5         
#>  Min.   :-16.069  
#>  1st Qu.: -8.519  
#>  Median : -8.113  
#>  Mean   : -8.241  
#>  3rd Qu.: -7.829  
#>  Max.   : -7.265  
#>        V1                V2                V3                 V4         
#>  Min.   :0.06622   Min.   :0.06211   Min.   :0.006447   Min.   :0.02927  
#>  1st Qu.:0.68746   1st Qu.:0.67060   1st Qu.:0.562012   1st Qu.:0.10470  
#>  Median :0.81054   Median :0.79601   Median :0.707682   Median :0.13142  
#>  Mean   :0.78220   Mean   :0.76862   Mean   :0.688361   Mean   :0.13573  
#>  3rd Qu.:0.90889   3rd Qu.:0.89728   3rd Qu.:0.836746   3rd Qu.:0.16257  
#>  Max.   :0.99995   Max.   :0.99998   Max.   :0.999951   Max.   :0.32660  
#>        V5         
#>  Min.   :-13.019  
#>  1st Qu.: -8.518  
#>  Median : -8.111  
#>  Mean   : -8.239  
#>  3rd Qu.: -7.829  
#>  Max.   : -7.274  
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
#>        V1                V2                V3                 V4         
#>  Min.   :0.03093   Min.   :0.01757   Min.   :0.001359   Min.   :0.01908  
#>  1st Qu.:0.68424   1st Qu.:0.67174   1st Qu.:0.560334   1st Qu.:0.10266  
#>  Median :0.80790   Median :0.79561   Median :0.705529   Median :0.12861  
#>  Mean   :0.78038   Mean   :0.76890   Mean   :0.686263   Mean   :0.13308  
#>  3rd Qu.:0.90730   3rd Qu.:0.89667   3rd Qu.:0.833137   3rd Qu.:0.15897  
#>  Max.   :0.99997   Max.   :0.99998   Max.   :0.999933   Max.   :0.36518  
#>        V5         
#>  Min.   :-12.690  
#>  1st Qu.: -8.522  
#>  Median : -8.113  
#>  Mean   : -8.238  
#>  3rd Qu.: -7.823  
#>  Max.   : -7.276  
#>        V1                V2               V3                 V4         
#>  Min.   :0.05163   Min.   :0.0469   Min.   :0.003935   Min.   :0.01756  
#>  1st Qu.:0.68257   1st Qu.:0.6737   1st Qu.:0.559975   1st Qu.:0.10127  
#>  Median :0.81103   Median :0.7985   Median :0.710753   Median :0.12984  
#>  Mean   :0.78060   Mean   :0.7697   Mean   :0.687697   Mean   :0.13294  
#>  3rd Qu.:0.90604   3rd Qu.:0.8994   3rd Qu.:0.837917   3rd Qu.:0.15818  
#>  Max.   :0.99993   Max.   :1.0000   Max.   :0.999967   Max.   :0.33020  
#>        V5         
#>  Min.   :-13.368  
#>  1st Qu.: -8.539  
#>  Median : -8.124  
#>  Mean   : -8.251  
#>  3rd Qu.: -7.836  
#>  Max.   : -7.257  
#>        V1                V2                 V3                  V4         
#>  Min.   :0.02559   Min.   :0.004825   Min.   :0.0004584   Min.   :0.02514  
#>  1st Qu.:0.68012   1st Qu.:0.667999   1st Qu.:0.5610183   1st Qu.:0.10480  
#>  Median :0.80770   Median :0.801639   Median :0.7095216   Median :0.12999  
#>  Mean   :0.77847   Mean   :0.770320   Mean   :0.6883072   Mean   :0.13435  
#>  3rd Qu.:0.90716   3rd Qu.:0.899039   3rd Qu.:0.8387168   3rd Qu.:0.16036  
#>  Max.   :0.99998   Max.   :0.999988   Max.   :0.9997415   Max.   :0.33959  
#>        V5         
#>  Min.   :-12.115  
#>  1st Qu.: -8.514  
#>  Median : -8.120  
#>  Mean   : -8.238  
#>  3rd Qu.: -7.831  
#>  Max.   : -7.270  
#>        V1                V2                 V3                V4         
#>  Min.   :0.02657   Min.   :0.001194   Min.   :0.01719   Min.   :0.01838  
#>  1st Qu.:0.68027   1st Qu.:0.673044   1st Qu.:0.56302   1st Qu.:0.10445  
#>  Median :0.80647   Median :0.797460   Median :0.70861   Median :0.13058  
#>  Mean   :0.77829   Mean   :0.769806   Mean   :0.68872   Mean   :0.13466  
#>  3rd Qu.:0.90772   3rd Qu.:0.897458   3rd Qu.:0.83349   3rd Qu.:0.16048  
#>  Max.   :0.99999   Max.   :0.999974   Max.   :0.99992   Max.   :0.33365  
#>        V5         
#>  Min.   :-12.848  
#>  1st Qu.: -8.510  
#>  Median : -8.118  
#>  Mean   : -8.240  
#>  3rd Qu.: -7.828  
#>  Max.   : -7.283  
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
#> [1,] 0.00000 0.27939 0.14831
#> [2,] 0.19691 0.00000 0.26440
#> [3,] 0.05346 0.16891 0.00000
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
