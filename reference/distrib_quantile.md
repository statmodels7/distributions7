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
