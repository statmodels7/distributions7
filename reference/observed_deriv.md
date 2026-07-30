# Observed Derivatives of a Given Order

Routes to
[`distrib_gradient`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.md),
[`distrib_hessian`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md),
[`distrib_deriv3`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
or
[`distrib_deriv4`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
according to `order`, so that code working at an order fixed only at run
time does not have to branch.

## Usage

``` r
observed_deriv(distrib, y, theta, order)
```

## Arguments

- distrib:

  An object inheriting from class `"distrib"`.

- y:

  A numeric vector of observations.

- theta:

  A named list of parameters.

- order:

  The derivative order, 1 to 4.

## Value

A named list of derivative component vectors, keyed as
[`hess_names`](https://statmodels7.github.io/distributions7/reference/hess_names.md)
at order 2 and
[`deriv_names`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
above it.
