# Pseudo-Huber Log-CDF Hessian

Closed form in the location-scale block; the three components touching
the shape are differenced. The method is
[`partial_loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_hess_cdf.md)
itself, shared with the Student t and the skew t.

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu`, `sigma` (positive) and `nu`
  (positive), each a numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

## Value

A named list of six numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `q` recycled against `theta`.

## Details

The saving is larger here than for the Student t, this family's
distribution function being a quadrature: the closed block reads the
density and its response derivative, where differencing it would run the
quadrature four times per component.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\\nu \> 0\\ the
shape, \\z = (q-\mu)/\sigma\\ and \\f\\ the density.

## See also

[`partial_loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/partial_loc_scale_hess_cdf.md)
for the shared body;
[`distrib_grad_cdf.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.PseudoHuberDistrib.md)
for the first order;
[`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

## Examples

``` r
d <- pseudohuber_distrib()
th <- list(mu = 0.3, sigma = 1.2, nu = 4)

distrib_hess_cdf(d, c(-1, 2), th)$mu_mu
#> [1] -0.1121827 -0.1037853
```
