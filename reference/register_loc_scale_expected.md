# Register the Location-Scale Expected Information on a Family

Registers
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
on a location-scale class, each computed by
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md),
and
[`expected_hessian_by_quadrature()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_by_quadrature.md)
answering `TRUE`.

## Usage

``` r
register_loc_scale_expected(cls)
```

## Arguments

- cls:

  An S7 class whose first two parameters are a location and a scale.

## Value

`NULL`, invisibly; called for its side effect.

## Details

The `approx` and `nsim` arguments of the three generics are accepted and
ignored, the expectation being taken by the rule of
[`loc_scale_rule()`](https://statmodels7.github.io/distributions7/reference/loc_scale_rule.md)
and not by one of the approximations `approx` selects between.
