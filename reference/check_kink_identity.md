# Stop at the Hessian Identity of a Family With a Kink

Signals an error when the Hessian identity is requested at order 2 for a
family that has a non-smooth parameter, because the observed second
derivative then lacks the point mass at the kink.

## Usage

``` r
check_kink_identity(distrib, order, who)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- order:

  The derivative order.

- who:

  The name of the constructor, used in the message.

## Value

`NULL`, invisibly.

## See also

[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
