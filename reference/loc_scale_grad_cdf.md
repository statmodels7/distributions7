# Location-Scale CDF Gradient

This function is the
[`distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.md)
body that the location-scale families share: it computes the factors of
[`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
at order 1 and puts them on the requested tail and scale by
[`loc_scale_cdf_route()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_route.md).
The Gaussian, the logistic, the Laplace and the Gumbel register this
function itself, so their pages document one piece of arithmetic.

## Usage

``` r
loc_scale_grad_cdf(distrib, q, theta, lower.tail = TRUE, log = TRUE, ...)
```

## Arguments

- distrib:

  An object inheriting from `distrib` whose first two parameters are a
  location and a scale.

- q:

  A numeric vector of quantiles.

- theta:

  A named list of parameters, in which the location comes first and the
  scale second.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, one per parameter, in which each
vector has the length of `q` recycled against `theta`.

## See also

[`loc_scale_cdf_deriv()`](https://statmodels7.github.io/distributions7/reference/loc_scale_cdf_deriv.md)
for the formulas;
[`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md)
for the second order;
[`distrib_grad_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Gaussian1Distrib.md)
for one of the four registrations.
