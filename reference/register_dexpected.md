# Register a Family's Derivatives of the Expected Information

Registers methods of
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
for one class, both reading one kernel on the parameter scale through
[`dexpected_analytic()`](https://statmodels7.github.io/distributions7/reference/dexpected_analytic.md).

## Usage

``` r
register_dexpected(cls, kernel)
```

## Arguments

- cls:

  An S7 class.

- kernel:

  A function of `(distrib, y, theta, order, threads)` returning the
  parameter-scale components keyed as
  [`dexpected_names()`](https://statmodels7.github.io/distributions7/reference/dexpected_names.md)
  at order 1 and
  [`d2expected_names()`](https://statmodels7.github.io/distributions7/reference/d2expected_names.md)
  at order 2.

## Value

`NULL`, invisibly; called for the registration.
