# Location-Scale CDF Hessian

The
[`distrib_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.md)
body the location-scale families share: the factors of
[`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
at both orders, put on the requested tail and scale by
[`loc_scale_cdf_route()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_route.md).
Both orders are needed even when only the Hessian was requested, the
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

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

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
