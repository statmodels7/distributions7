# The Order a Fifth Derivative Is Differenced From

Returns 4 where
[`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md)
is `TRUE`; 0 for a wrapper whose parent has no analytic fourth order, a
wrapper being a class that overrides
[`has_exact_deriv4()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv4.md)
to read its parent; and otherwise the highest order the family registers
itself, as
[`analytic_order()`](https://statmodels7.github.io/distributions7/reference/analytic_order.md)
gives it up to 3.

## Usage

``` r
deriv5_base(x)
```

## Arguments

- x:

  A distribution object.

## Value

An integer from 0 to 4.

## See also

[`numerical_deriv5()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv5.md),
which passes it to
[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md).
