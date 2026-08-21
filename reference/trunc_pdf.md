# Density of a Truncated Distribution

The parent's density divided by \\Z\\, and zero outside the interval.

## Usage

``` r
trunc_pdf(distrib, y, theta, log = FALSE, ...)
```

## Arguments

- distrib:

  A truncated distribution object.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- log:

  Logical; if `TRUE`, returns the log-density.

## Value

A numeric vector.

## Details

One of the shared method bodies. Truncation treats the two kinds of
parent identically once
[`trunc_constants`](https://statmodels7.github.io/distributions7/reference/trunc_constants.md)
has resolved the one place they differ, so these bodies are written once
and registered on both classes.
