# Quantile Function

Evaluates the quantile function for a given distribution.

## Usage

``` r
distrib_quantile(distrib, p, theta, ...)
```

## Arguments

- distrib:

  A distribution object inheriting from the `distrib` class.

- p:

  A numeric vector of probabilities.

- theta:

  A named list (or named numeric vector) of distribution parameters. If
  unnamed, parameters are taken in the order of `distrib@params`.

- ...:

  Additional arguments passed to the specific method (e.g.,
  `lower.tail`, `log.p`).

## Value

A numeric vector of quantiles.

## Examples

``` r
distrib_quantile(gaussian1_distrib(), c(0.025, 0.5, 0.975), list(mu = 0, sigma = 1))
#> [1] -1.959964  0.000000  1.959964
distrib_quantile(poisson_distrib(), c(0.1, 0.9), list(mu = 2))
#> [1] 0 4
```
