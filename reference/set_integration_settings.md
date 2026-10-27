# Set the Integration Settings of an Approximation

Makes
[`quad_rows()`](https://statmodels7.github.io/distributions7/reference/quad_rows.md)
and
[`series_rows()`](https://statmodels7.github.io/distributions7/reference/series_rows.md)
read their tolerances and limits from `method`, and returns the settings
that were in force so that the caller can restore them with
[`on.exit()`](https://rdrr.io/r/base/on.exit.html).

## Usage

``` r
set_integration_settings(method)
```

## Arguments

- method:

  An
  [ExpectIntegrate](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
  object, or `NULL` for the defaults.

## Value

The previous settings, invisibly: an
[ExpectIntegrate](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md)
object or `NULL`.

## See also

[`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md),
[`integration_settings()`](https://statmodels7.github.io/distributions7/reference/integration_settings.md)
