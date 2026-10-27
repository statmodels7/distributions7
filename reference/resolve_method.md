# Resolve the Default Method of an Expected Derivative

Returns `method` when it is an object and, when it is `NULL`, the
default of the generics
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the order:
[`expect_opg()`](https://statmodels7.github.io/distributions7/reference/expect_opg.md)
at order 2 and
[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md)
with `identity = "hessian"` above it.

## Usage

``` r
resolve_method(method, order)
```

## Arguments

- method:

  An object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
  or `NULL`.

- order:

  The derivative order, 2 to 4.

## Value

An object inheriting from
[ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).

## See also

[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
