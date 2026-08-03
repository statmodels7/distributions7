# Refuse a Quantity That Has No Multivariate Counterpart

Raises the error a multivariate distribution gives for the
one-dimensional quantities: the distribution function, the quantile
function and the expectation by quadrature.

## Usage

``` r
mv_refuse(distrib, what, why)
```

## Arguments

- distrib:

  A
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md)
  object.

- what:

  The name of the quantity.

- why:

  One sentence saying what is missing.

## Value

Never returns; raises an error.
