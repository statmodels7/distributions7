# The Closed-Form Higher Derivatives of a Multivariate Gaussian

Shared engine for the third and fourth derivatives: enumerates the
parameter tuples the way
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
does, splits each into mean and structure indices, and reads the
surviving cases off the gaussian's algebra.

## Usage

``` r
mvg_higher(distrib, y, theta, order)
```

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- y:

  An \\n imes p\\ matrix of observations.

- theta:

  A named list of parameters.

- order:

  3 or 4.

## Value

A named list of derivative component vectors.
