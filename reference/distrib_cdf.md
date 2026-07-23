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
