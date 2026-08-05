# Precision Derivative Tensors of a Multivariate Gaussian

The precision's derivative tensors in the matrix parameter's free
values, orders 1 to 4, keyed by index tuple. For a precision
parametrization they are the parameter's own derivatives; for a
covariance they follow from repeated differentiation of the inverse, so
no expanded formula is transcribed and no term can be dropped.

## Usage

``` r
mvg_ptensors(distrib, theta, order)
```

## Arguments

- distrib:

  A
  [`MvGaussianDistrib`](https://statmodels7.github.io/distributions7/reference/MvGaussianDistrib.md)
  object.

- theta:

  A named list of parameters.

- order:

  The highest order wanted.

## Value

A list with the accessor `get`, the log-determinant sign and the pieces.
