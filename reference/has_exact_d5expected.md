# Does a Family Supply Its Own Fifth Derivative of the Information?

Returns `TRUE` where the
[`distrib_d5expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d5expected_hessian.md)
method that dispatches is registered on a class other than the package's
base classes, whose method is the stencil of
[`numerical_d5expected()`](https://statmodels7.github.io/distributions7/reference/numerical_d5expected.md).

## Usage

``` r
has_exact_d5expected(x)
```

## Arguments

- x:

  A distribution object.

## Value

A single logical.

## See also

[`distrib_d5expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d5expected_hessian.md),
which selects its route on the link scale with it.
