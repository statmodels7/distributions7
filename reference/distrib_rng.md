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
