# The Order from Which a Sixth Derivative Is Differenced

Returns 5 where the family's fifth derivative is analytic, as
[`has_exact_deriv5()`](https://statmodels7.github.io/distributions7/reference/has_exact_deriv5.md)
reports it, and otherwise the order that
[`deriv5_base()`](https://statmodels7.github.io/distributions7/reference/deriv5_base.md)
returns.

## Usage

``` r
deriv6_base(x)
```

## Arguments

- x:

  A distribution object.

## Value

An integer from 0 to 5.

## See also

[`numerical_deriv6()`](https://statmodels7.github.io/distributions7/reference/numerical_deriv6.md),
which passes it to
[`tensor_derivatives()`](https://statmodels7.github.io/distributions7/reference/tensor_derivatives.md).
