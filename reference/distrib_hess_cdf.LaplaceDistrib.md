# Laplace Log-CDF Hessian

Closed form, from the same location-scale structure, through
[`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md).
The response derivative it reads is \\\ell_y =
-\mathrm{sign}(z)/\sigma\\, which jumps at \\q = \mu\\, so the second
derivatives of \\F\\ jump there too. That is a property of the law and
not a defect: the Laplace is the toolkit's non-regular family, and its
location has no second derivative at the kink.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `sigma` (positive), each a
  numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

## Value

A named list of three numeric vectors keyed as
[`hess_names()`](https://statmodels7.github.io/distributions7/reference/hess_names.md),
each the length of `q` recycled against `theta`.

## Details

Away from \\q = \mu\\ the three components are the ordinary
location-scale ones. At \\q = \mu\\ exactly, `sign(0)` is 0 and the
returned `mu_mu` is the average of the two one-sided limits, which is
the value the sign convention gives and is reported as it stands.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\\ell_y = \partial\log
f/\partial y\\.

## See also

[`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md)
for the shared body;
[`distrib_grad_cdf.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.LaplaceDistrib.md),
which is continuous at the kink;
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

## Examples

``` r
d <- laplace_distrib()
th <- list(mu = 0.3, sigma = 1.2)

# The second derivative in the location jumps across q = mu.
distrib_hess_cdf(d, 0.3 + c(-1e-6, 1e-6), th, log = FALSE)$mu_mu
#> [1]  0.3472219 -0.3472219

# Away from the kink it agrees with a central difference of the cdf.
q <- c(-1, 2)
exact <- distrib_hess_cdf(d, q, th, log = FALSE)
fd <- numerical_cdf_deriv(d, q, th, order = 2)
max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
#> [1] 1.064456e-07
```
