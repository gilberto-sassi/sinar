
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sinar

<!-- badges: start -->

<!-- badges: end -->

The goal of sinar is to implement the Conditional Least Square method
for the Spatial non-negative Integer-valued Autoregressive
\(SINAR(1,1)\).

## Installation

And the development version from
[GitHub](https://github.com/gilberto-sassi/) with:

``` r
# install.packages("devtools")
devtools::install_github("gilberto-sassi/sinar")
```

## Example: simulated case

``` r
library(sinar)

## Simulated data matrix from SINAR(1,1) with Poison(5) innovation
matrix_simulated <- sinar_pois(15, 15, 0.2, 0.2, 0.4, 5)

## Conditional Least Square (CLS) estimates
cls(matrix_simulated)
#>       a10       a01       a11        mu 
#> 0.2658153 0.1906724 0.3980458 3.9129349

## Covariance matrix of CLS estimates
emp_cov(matrix_simulated)
#>               a10          a01           a11          mu
#> a10  0.0034154179  0.000228933 -0.0009146087 -0.06524879
#> a01  0.0002289330  0.003859600 -0.0016068704 -0.05990849
#> a11 -0.0009146087 -0.001606870  0.0031029738 -0.01236036
#> mu  -0.0652487913 -0.059908492 -0.0123603642  3.35664781
```

## Example: real dataset (nematodoes)

``` r
library(sinar)

## Nematodes counting datasets
data("nematodes")

## Conditional Least Square (CLS) estimates
cls(nematodes)
#>        a10        a01        a11         mu 
#> 0.20664577 0.33147378 0.04523086 2.14476453

## Covariance matrix of CLS estimates
emp_cov(nematodes)
#>               a10           a01          a11           mu
#> a10  0.0111169222 -0.0009999304 -0.003310576 -0.017278481
#> a01 -0.0009999304  0.0082946407 -0.001503724 -0.009838536
#> a11 -0.0033105760 -0.0015037242  0.004507501  0.004049939
#> mu  -0.0172784806 -0.0098385364  0.004049939  0.268045835
```

## Example: real dataset (carabidae)

``` r
library(sinar)

## Carabidae counting dataset
data("carabidae")

## Conditional Least Square (CLS) estimates
cls(carabidae)
#>        a10        a01        a11         mu 
#> 0.14595392 0.12725313 0.08798513 9.10361759

## Covariance matrix of CLS estimates
emp_cov(carabidae)
#>              a10          a01          a11          mu
#> a10  0.014484776 -0.003141815 -0.005525906 -0.06795645
#> a01 -0.003141815  0.014365625 -0.001265544 -0.11558802
#> a11 -0.005525906 -0.001265544  0.023795735 -0.25417404
#> mu  -0.067956449 -0.115588024 -0.254174036  7.22525572
```
