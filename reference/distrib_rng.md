# Random Number Generator

Generates random variates from the given distribution.

## Usage

``` r
distrib_rng(distrib, n, theta, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- n:

  Number of observations to generate.

- theta:

  A named list (or named numeric vector) of distribution parameters. If
  unnamed, parameters are taken in the order of `distrib@params`.

- ...:

  Additional arguments passed to the specific method.

## Value

A numeric vector of `n` draws for a univariate distribution, and an \\n
\times p\\ matrix for a multivariate one.

## Examples

``` r
set.seed(1)
distrib_rng(gaussian1_distrib(), 5, list(mu = 0, sigma = 1))
#> [1] -0.6264538  0.1836433 -0.8356286  1.5952808  0.3295078
distrib_rng(mvgaussian_distrib(2), 3, list(mu1 = 0, mu2 = 0,
  sigma_log_L1 = 0, sigma_log_L2 = 0, sigma_L2.1 = 0.5))
#>              v1          v2
#> [1,] -0.8204684  0.16554716
#> [2,]  0.4874291 -0.06167386
#> [3,]  0.7383247  1.88094352
```
