# Location-Scale CDF Hessian

The
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
body the location-scale families share:
[`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
at both orders, put on the requested tail and scale by
[`cdf_tail_scale()`](https://statmodels7.github.io/distributions7/reference/cdf_tail_scale.md).
Both orders are needed even when only the Hessian was asked for, the
log-scale correction reading the first derivatives.

## Usage

``` r
loc_scale_hess_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...)
```

## Arguments

- distrib:

  An object inheriting from `distrib` whose first two parameters are a
  location and a scale.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters, the location first and the scale second.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `q` recycled against `theta`. The gradient is not
returned alongside.

## See also

[`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
for the formulas;
[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
for the first order;
[`distrib_hess_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.Gaussian1Distrib.md)
for one of the four registrations.
