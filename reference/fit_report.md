# Standard Errors and Intervals of a Fit for a Chosen Information

Returns the quantities of
[`fit_inference()`](https://statmodels7.github.io/distributions7/reference/fit_inference.md)
for the information selected by an
[Information](https://statmodels7.github.io/distributions7/reference/Information.md)
object, or by the default of the family when `information` is `NULL`.
The fit stores the quantities of the default information
([`default_information()`](https://statmodels7.github.io/distributions7/reference/default_information.md)),
and they are returned when the selected information is the default;
otherwise the link-scale variance matrix is recomputed by
[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md)
and carried through
[`fit_inference()`](https://statmodels7.github.io/distributions7/reference/fit_inference.md).

## Usage

``` r
fit_report(object, information = NULL)
```

## Arguments

- object:

  A
  [`distrib_fit()`](https://statmodels7.github.io/distributions7/reference/distrib_fit_class.md)
  object.

- information:

  `NULL`, an
  [`observed_info()`](https://statmodels7.github.io/distributions7/reference/observed_info.md)
  or an
  [`expected_info()`](https://statmodels7.github.io/distributions7/reference/expected_info.md)
  object.

## Value

A list as returned by
[`fit_inference()`](https://statmodels7.github.io/distributions7/reference/fit_inference.md).

## See also

[`vcov.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/vcov.distrib_fit.md),
[`print.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/print.distrib_fit.md),
[`confint.distrib_fit()`](https://statmodels7.github.io/distributions7/reference/confint.distrib_fit.md),
[`mv_summary()`](https://statmodels7.github.io/distributions7/reference/mv_summary.md)
