# Does a Family Supply Its Own Sixth Derivative?

Returns `TRUE` where the
[`distrib_deriv6()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv6.md)
method that dispatches is registered on a class other than the package's
base classes, whose method is the stencil of
[`numerical_deriv6()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv6.md).
A wrapper registers no sixth-order method, so the result is `FALSE` for
it.

## Usage

``` r
has_exact_deriv6(x)
```

## Arguments

- x:

  A distribution object.

## Value

A single logical.

## See also

[`distrib_deriv6()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv6.md),
which selects its route on the link scale with it.
