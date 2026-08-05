# Every Component of a Beta-Binomial Derivative

Assembles the components of a derivative of the given order, named as
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them.

## Usage

``` r
betabinom2_derivs(y, a, b, n, order, params)
```

## Arguments

- y:

  A numeric vector of observations.

- a, b:

  The two shapes.

- n:

  The size.

- order:

  The derivative order.

- params:

  The parameter names.

## Value

A named list of component vectors.

## See also

[`betabinom2_distrib`](https://statmodels7.github.io/distributions7/reference/betabinom2_distrib.md)
