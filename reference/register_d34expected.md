# Register a Family's Third and Fourth Derivatives of the Expected Information

Registers methods of
[`distrib_d3expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d3expected_hessian.md)
and
[`distrib_d4expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d4expected_hessian.md)
for one class, both reading one kernel on the parameter scale through
[`d34expected_analytic()`](https://statmodels7.github.io/distributions7/reference/d34expected_analytic.md).

## Usage

``` r
register_d34expected(cls, kernel)
```

## Arguments

- cls:

  An S7 class.

- kernel:

  A function of `(distrib, y, theta, order, threads)` returning, for
  `order` 3 or 4, the parameter-scale components in the order of
  [`d3expected_names()`](https://statmodels7.github.io/distributions7/reference/d3expected_names.md)
  or
  [`d4expected_names()`](https://statmodels7.github.io/distributions7/reference/d4expected_names.md);
  names are set by the caller.

## Value

`NULL`, invisibly; called for the registration.
