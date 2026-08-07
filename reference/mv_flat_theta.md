# Require Scalar Parameters

Rejects a `theta` whose components are not single numbers.

## Usage

``` r
mv_flat_theta(distrib, theta)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- theta:

  A named list of parameters.

## Value

A numeric vector of the parameter values, in declaration order.

## Details

The multivariate families of this package take one parameter value for
the whole sample. Vectorized parameters are what a regression supplies,
and a distribution that accepted them would be doing the model layer's
work with none of its bookkeeping.
