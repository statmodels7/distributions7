# Register the Location-Scale Expected Information on a Family

Registers
[`distrib_expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.md),
[`distrib_dexpected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_dexpected_hessian.md)
and
[`distrib_d2expected_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_d2expected_hessian.md)
on a location-scale class, each of which is computed by
[`loc_scale_expected()`](https://statmodels7.github.io/distributions7/reference/loc_scale_expected.md),
and registers
[`expected_hessian_costly()`](https://statmodels7.github.io/distributions7/reference/expected_hessian_costly.md)
so that it returns `TRUE`.

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

The `method` arguments of the three methods receive `NULL`, because the
generic applies an object before dispatching, and the expectation is
then taken by the rule of
[`loc_scale_rule()`](https://statmodels7.github.io/distributions7/reference/loc_scale_rule.md).
