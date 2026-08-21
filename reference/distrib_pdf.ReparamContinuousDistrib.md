# Density of a Reparametrized Distribution

The parent's density at the mapped parameters.

## Usage

``` r
reparam_pdf(distrib, y, theta, log = FALSE, ...)
```

## Arguments

- distrib:

  A reparametrized distribution.

- y:

  A numeric vector of observations.

- theta:

  A named list of the new parameters.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector.

## See also

[`reparametrize`](https://statmodels7.github.io/distributions7/reference/reparametrize.md)
