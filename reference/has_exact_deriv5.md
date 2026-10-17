# Is a Family's Fifth Derivative Analytic?

Returns `TRUE` where the
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
method that dispatches is registered on a class other than the base
class, whose method is the stencil of
[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md).
For a wrapper registered in this file, the result is that of the family
that the wrapper wraps, since the wrapper's fifth derivative is analytic
exactly when that family's is.

## Usage

``` r
has_exact_deriv5(x)
```

## Arguments

- x:

  A distribution object.

## Value

A single logical.
