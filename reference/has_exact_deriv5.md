# Is a Family's Fifth Derivative Analytic?

`TRUE` where the
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
method that dispatches is registered on a class other than the base
class, whose method is the stencil of
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md).
A wrapper registered in this file answers for the family it wraps, since
its fifth derivative is analytic exactly when that family's is.

## Usage

``` r
has_exact_deriv5(x)
```

## Arguments

- x:

  A distribution object.

## Value

A single logical.
