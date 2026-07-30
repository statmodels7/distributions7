# Random Generation From a Truncated Distribution

Draws by inverse transform through
[`trunc_quantile`](https://statmodels7.github.io/distributions7/reference/trunc_quantile.md),
which is exact and needs no rejection.

## Usage

``` r
trunc_rng(distrib, n, theta)
```

## Arguments

- distrib:

  A truncated distribution object.

- n:

  The number of draws.

- theta:

  A named list of parameters.

## Value

A numeric vector of draws.
