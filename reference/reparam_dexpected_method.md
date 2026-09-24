# A Method Body for the Reparametrized Families

Returns the method of
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
(order 1) or
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
(order 2) registered on the reparametrized classes, reading
[`reparam_dexpected()`](https://statmodels7.github.io/distributions7/reference/reparam_dexpected.md).

## Usage

``` r
reparam_dexpected_method(order)
```

## Arguments

- order:

  `1L` or `2L`.

## Value

A function with the generics' signature.
