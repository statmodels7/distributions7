# Cumulative Distribution Function

Evaluates the cumulative distribution function (CDF) for a given
distribution.

## Usage

``` r
distrib_cdf(distrib, q, theta, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- q:

  A numeric vector of quantiles.

- theta:

  A named list (or named numeric vector) of distribution parameters. If
  unnamed, parameters are taken in the order of `distrib@params`.

- ...:

  Additional arguments passed to the specific method (e.g.,
  `lower.tail`, `log.p`).

## Value

A numeric vector of cumulative probabilities.

## Examples

``` r
distrib_cdf(gaussian_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> [1] 0.1586553 0.5000000 0.8413447
distrib_cdf(poisson_distrib(), 0:3, list(mu = 2), lower.tail = FALSE)
#> [1] 0.8646647 0.5939942 0.3233236 0.1428765
```
