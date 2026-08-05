# The First p Parameters, Read as a Location

The helper the elliptical families implement
[`mv_location`](https://statmodels7.github.io/distributions7/reference/mv_location.md)
with: the first \\p\\ entries of the flat parameter vector, labeled by
coordinate.

## Usage

``` r
mv_leading_location(distrib, theta)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters, already aligned.

## Value

A named numeric vector of length \\p\\.
