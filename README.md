
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Stratified Sampling

Integrating a stratified structure in the population in a sampling
design can considerably reduce the variance of the Horvitz-Thompson
estimator. We propose in this package different methods to handle the
selection of a balanced sample in stratified population. For more
details see Raphaël Jauslin, Esther Eustache and Yves Tillé (2021)
<https://arxiv.org/abs/2101.05568>.

## Build

[![Build
Status](https://travis-ci.org/RJauslin/StratifiedSampling.svg?branch=master)](https://travis-ci.org/RJauslin/StratifiedSampling)

## Installation

<!-- ### CRAN version -->

<!-- ``` -->

<!-- install.packages("WaveSampling") -->

<!-- ``` -->

### Latest version

You can install the latest version of the package `StratifiedSampling`
with the following command:

``` r
# install.packages("devtools")
devtools::install_github("Rjauslin/StratifiedSampling")
```

## Simple example

This basic example shows you how to set up a stratified sampling design.
The example is done on the `swissmunicipalities` dataset from the
package `sampling`.

``` r
library(sampling)
library(StratifiedSampling)
#> Loading required package: Matrix
#> 
#> Attaching package: 'StratifiedSampling'
#> The following object is masked from 'package:base':
#> 
#>     choose

data(swissmunicipalities)
swiss <- swissmunicipalities
X <- cbind(swiss$HApoly,
        swiss$Surfacesbois,
        swiss$P00BMTOT,
        swiss$P00BWTOT,
        swiss$POPTOT,
        swiss$Pop020,
        swiss$Pop2040,
        swiss$Pop4065,
        swiss$Pop65P,
        swiss$H00PTOT )

X <- X[order(swiss$REG),]
strata <- swiss$REG[order(swiss$REG)]
```

Strata are NUTS region of the Switzerland. Inclusion probabilities `pik`
is set up equal within strata and such that the sum of the inclusion
probabilities within strata is equal to 80.

``` r
pik <- sampling::inclusionprobastrata(strata,rep(80,7))
```

It remains to use the function `stratifiedcube()`.

``` r
s <- stratifiedcube(X,strata,pik)
```

We can check that we have correctly selected the sample. It is balanced
and have the right number of units selected in each stratum.

``` r
head(s)
#> [1] 0 0 0 0 0 0

sum(s)
#> [1] 560
t(X/pik)%*%s
#>          [,1]
#>  [1,] 4035984
#>  [2,] 1256920
#>  [3,] 3439487
#>  [4,] 3575885
#>  [5,] 7015372
#>  [6,] 1615537
#>  [7,] 2049134
#>  [8,] 2281516
#>  [9,] 1069184
#> [10,] 2974870
t(X/pik)%*%pik
#>          [,1]
#>  [1,] 3998831
#>  [2,] 1270996
#>  [3,] 3567567
#>  [4,] 3720443
#>  [5,] 7288010
#>  [6,] 1665613
#>  [7,] 2141059
#>  [8,] 2362332
#>  [9,] 1119006
#> [10,] 3115399

Xcat <- disj(strata)

t(Xcat)%*%s
#>      [,1]
#> [1,]   80
#> [2,]   80
#> [3,]   80
#> [4,]   80
#> [5,]   80
#> [6,]   80
#> [7,]   80
t(Xcat)%*%pik
#>      [,1]
#> [1,]   80
#> [2,]   80
#> [3,]   80
#> [4,]   80
#> [5,]   80
#> [6,]   80
#> [7,]   80
```
