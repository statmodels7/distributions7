# Cauchy Log-CDF Hessian

Computes the Hessian in closed form, from the same location-scale
structure, through
[`cauchy_cdf_method()`](https://statmodels7.github.io/distributions7/reference/cauchy_cdf_method.md),
which computes each component as one term in \\Y^{-2}\\, \\Y = q - \mu +
\mathrm{i}\sigma\\. With \\\ell_y = \partial\log f/\partial y\\, which
for a Cauchy is \\-2z/\\\sigma(1+z^2)\\\\, \$\$\frac{\partial^2
F}{\partial\mu^2} = f\\\ell_y, \qquad \frac{\partial^2
F}{\partial\mu\\\partial\sigma} = f\left(z\ell_y +
\frac{1}{\sigma}\right), \qquad \frac{\partial^2 F}{\partial\sigma^2} =
f\left(z^2\ell_y + \frac{2z}{\sigma}\right).\$\$ As with the gradient,
these are exact although the family has no moments.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `sigma` (positive); each is a
  numeric vector of length 1 or `n`.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability.

## Value

A named list of three numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
in which each vector has the length of `q` recycled against `theta`.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\\ell_y = \partial\log
f/\partial y\\.

## See also

[`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md)
for the shared body;
[`distrib_grad_cdf.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.CauchyDistrib.md)
for the first order;
[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

## Examples

``` r
d <- cauchy_distrib()
th <- list(mu = 0.3, sigma = 1.2)
q <- c(-1, 0.5, 2)

# The result is compared with a central difference of the cdf, which shares
# no arithmetic with it.
exact <- distrib_hess_cdf(d, q, th, log = FALSE)
fd <- numerical_cdf_deriv(d, q, th, order = 2)
max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
#> [1] 8.886859e-08
```
