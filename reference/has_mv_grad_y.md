# Whether a Multivariate Family Implements Its Response Derivatives

`TRUE` when
[`distrib_grad_y`](https://statmodels7.github.io/distributions7/reference/distrib_grad_y.md)
comes from the family rather than from the base class, whose method
rejects.

## Usage

``` r
has_mv_grad_y(x)
```

## Arguments

- x:

  An object inheriting from class
  [`multivariate_distrib`](https://statmodels7.github.io/distributions7/reference/multivariate_distrib.md).

## Value

A single logical.

## See also

[`has_mv_support`](https://statmodels7.github.io/distributions7/reference/has_mv_support.md),
[`check_distrib`](https://statmodels7.github.io/distributions7/reference/check_distrib.md)
