# Whether an Expectation Method Replaces the Family's Own Route

Returns `TRUE` when `method` is an object. The generics
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
then compute the expectation with
[`route_expected()`](https://statmodels7.github.io/distributions7/reference/route_expected.md),
bypassing the method registered for the family, so that an explicit
method is used also by a family that writes the expectation in closed
form.

## Usage

``` r
explicit_method(distrib, method)
```

## Arguments

- distrib:

  An object inheriting from `distrib`.

- method:

  An object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md),
  or `NULL`.

## Value

A single logical.

## See also

[expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md)
