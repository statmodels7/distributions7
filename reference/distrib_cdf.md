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

## Details

\$\$F(q; \theta) = P(Y \le q),\$\$

the integral of the density up to \\q\\ for a continuous family and the
sum of the mass over the support points at or below \\q\\ for a discrete
one. `lower.tail = FALSE` returns \\1 - F(q; \theta)\\ and
`log.p = TRUE` its logarithm, both computed on the log scale where a
family provides one. Without a method the value comes from quadrature of
the density.

## Examples

``` r
distrib_cdf(gaussian1_distrib(), c(-1, 0, 1), list(mu = 0, sigma = 1))
#> [1] 0.1586553 0.5000000 0.8413447
distrib_cdf(poisson_distrib(), 0:3, list(mu = 2), lower.tail = FALSE)
#> [1] 0.8646647 0.5939942 0.3233236 0.1428765
```
