
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
devtools::install_github("rockwillck/CBNpath")
```

## Example

`CBNpath` provides an interface to CT-CBN, H-CBN, and B-CBN:

``` r
library(CBNpath)
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
#>        V1               V2                  V3                  V4         
#>  Min.   :0.1460   Min.   :0.0001472   Min.   :0.0000451   Min.   :0.01859  
#>  1st Qu.:0.8003   1st Qu.:0.2796521   1st Qu.:0.1418145   1st Qu.:0.11056  
#>  Median :0.8944   Median :0.4256849   Median :0.2839890   Median :0.13839  
#>  Mean   :0.8590   Mean   :0.4364513   Mean   :0.3334257   Mean   :0.14093  
#>  3rd Qu.:0.9548   3rd Qu.:0.5789421   3rd Qu.:0.4811127   3rd Qu.:0.16753  
#>  Max.   :1.0000   Max.   :0.9999941   Max.   :0.9998532   Max.   :0.32130  
#>        V5         
#>  Min.   :-12.872  
#>  1st Qu.: -8.440  
#>  Median : -8.028  
#>  Mean   : -8.146  
#>  3rd Qu.: -7.721  
#>  Max.   : -7.239  
#>        V1               V2                  V3                  V4         
#>  Min.   :0.1532   Min.   :0.0001293   Min.   :0.0000524   Min.   :0.02145  
#>  1st Qu.:0.7983   1st Qu.:0.2783313   1st Qu.:0.1413772   1st Qu.:0.10898  
#>  Median :0.8922   Median :0.4191844   Median :0.2789495   Median :0.13486  
#>  Mean   :0.8582   Mean   :0.4345782   Mean   :0.3323606   Mean   :0.13954  
#>  3rd Qu.:0.9537   3rd Qu.:0.5785409   3rd Qu.:0.4779362   3rd Qu.:0.16707  
#>  Max.   :1.0000   Max.   :0.9989193   Max.   :0.9999911   Max.   :0.34484  
#>        V5         
#>  Min.   :-13.118  
#>  1st Qu.: -8.451  
#>  Median : -8.048  
#>  Mean   : -8.162  
#>  3rd Qu.: -7.737  
#>  Max.   : -7.244  
#>        V1               V2                  V3                  V4         
#>  Min.   :0.1052   Min.   :0.0001035   Min.   :0.0000776   Min.   :0.02318  
#>  1st Qu.:0.8000   1st Qu.:0.2761990   1st Qu.:0.1404796   1st Qu.:0.10928  
#>  Median :0.8977   Median :0.4203094   Median :0.2843305   Median :0.13604  
#>  Mean   :0.8617   Mean   :0.4340512   Mean   :0.3330388   Mean   :0.13881  
#>  3rd Qu.:0.9547   3rd Qu.:0.5769578   3rd Qu.:0.4775278   3rd Qu.:0.16429  
#>  Max.   :1.0000   Max.   :0.9997298   Max.   :0.9995603   Max.   :0.38932  
#>        V5         
#>  Min.   :-12.863  
#>  1st Qu.: -8.438  
#>  Median : -8.027  
#>  Mean   : -8.145  
#>  3rd Qu.: -7.719  
#>  Max.   : -7.252  
#>        V1                V2                  V3                  V4         
#>  Min.   :0.02904   Min.   :0.0000646   Min.   :0.0000999   Min.   :0.02482  
#>  1st Qu.:0.79975   1st Qu.:0.2736313   1st Qu.:0.1406076   1st Qu.:0.10819  
#>  Median :0.89715   Median :0.4232942   Median :0.2785298   Median :0.13557  
#>  Mean   :0.85960   Mean   :0.4362567   Mean   :0.3292357   Mean   :0.13939  
#>  3rd Qu.:0.95584   3rd Qu.:0.5824931   3rd Qu.:0.4696677   3rd Qu.:0.16695  
#>  Max.   :1.00000   Max.   :0.9993357   Max.   :0.9994623   Max.   :0.34937  
#>        V5         
#>  Min.   :-13.040  
#>  1st Qu.: -8.450  
#>  Median : -8.035  
#>  Mean   : -8.157  
#>  3rd Qu.: -7.727  
#>  Max.   : -7.237  
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
#>        V1               V2                  V3                  V4         
#>  Min.   :0.1871   Min.   :0.0003938   Min.   :0.0000255   Min.   :0.02774  
#>  1st Qu.:0.8045   1st Qu.:0.2770992   1st Qu.:0.1395726   1st Qu.:0.10836  
#>  Median :0.8972   Median :0.4233658   Median :0.2787450   Median :0.13564  
#>  Mean   :0.8632   Mean   :0.4356416   Mean   :0.3277964   Mean   :0.13978  
#>  3rd Qu.:0.9574   3rd Qu.:0.5806498   3rd Qu.:0.4730893   3rd Qu.:0.16647  
#>  Max.   :1.0000   Max.   :0.9999167   Max.   :0.9999967   Max.   :0.34668  
#>        V5         
#>  Min.   :-12.057  
#>  1st Qu.: -8.423  
#>  Median : -8.019  
#>  Mean   : -8.142  
#>  3rd Qu.: -7.725  
#>  Max.   : -7.241  
#>        V1               V2                  V3                  V4        
#>  Min.   :0.1786   Min.   :0.0000179   Min.   :0.0000041   Min.   :0.0201  
#>  1st Qu.:0.7986   1st Qu.:0.2763707   1st Qu.:0.1425656   1st Qu.:0.1091  
#>  Median :0.8938   Median :0.4242553   Median :0.2809401   Median :0.1368  
#>  Mean   :0.8590   Mean   :0.4353905   Mean   :0.3325938   Mean   :0.1395  
#>  3rd Qu.:0.9558   3rd Qu.:0.5792342   3rd Qu.:0.4820772   3rd Qu.:0.1662  
#>  Max.   :1.0000   Max.   :0.9977905   Max.   :0.9999921   Max.   :0.3395  
#>        V5         
#>  Min.   :-12.501  
#>  1st Qu.: -8.459  
#>  Median : -8.034  
#>  Mean   : -8.156  
#>  3rd Qu.: -7.730  
#>  Max.   : -7.241  
#>        V1               V2                  V3                  V4         
#>  Min.   :0.1825   Min.   :0.0001759   Min.   :0.0000481   Min.   :0.03735  
#>  1st Qu.:0.8006   1st Qu.:0.2777662   1st Qu.:0.1411597   1st Qu.:0.11001  
#>  Median :0.8965   Median :0.4237951   Median :0.2794866   Median :0.13662  
#>  Mean   :0.8619   Mean   :0.4354211   Mean   :0.3310284   Mean   :0.13960  
#>  3rd Qu.:0.9580   3rd Qu.:0.5782369   3rd Qu.:0.4786070   3rd Qu.:0.16490  
#>  Max.   :1.0000   Max.   :0.9999084   Max.   :0.9996902   Max.   :0.33790  
#>        V5         
#>  Min.   :-12.469  
#>  1st Qu.: -8.430  
#>  Median : -8.019  
#>  Mean   : -8.141  
#>  3rd Qu.: -7.724  
#>  Max.   : -7.236  
#>        V1               V2                  V3                  V4         
#>  Min.   :0.1724   Min.   :0.0000759   Min.   :0.0000313   Min.   :0.02314  
#>  1st Qu.:0.7997   1st Qu.:0.2750210   1st Qu.:0.1421681   1st Qu.:0.10771  
#>  Median :0.8965   Median :0.4228841   Median :0.2835508   Median :0.13582  
#>  Mean   :0.8621   Mean   :0.4346064   Mean   :0.3314963   Mean   :0.13926  
#>  3rd Qu.:0.9550   3rd Qu.:0.5804758   3rd Qu.:0.4765703   3rd Qu.:0.16623  
#>  Max.   :1.0000   Max.   :0.9998170   Max.   :0.9999930   Max.   :0.32601  
#>        V5         
#>  Min.   :-12.704  
#>  1st Qu.: -8.435  
#>  Median : -8.027  
#>  Mean   : -8.143  
#>  3rd Qu.: -7.725  
#>  Max.   : -7.233  
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
#>         [,1]   [,2]    [,3]
#> [1,] 0.00000 0.2351 0.44870
#> [2,] 0.00150 0.0000 0.53604
#> [3,] 0.00011 0.0930 0.00000
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
