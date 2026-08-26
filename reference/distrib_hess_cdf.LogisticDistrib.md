# Logistic Log-CDF Hessian

Closed form, from the same location-scale structure, through
[`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md).
With \\\ell_y = \partial\log f/\partial y\\, which for a logistic is
\\(1-2F)/\sigma\\, \$\$\frac{\partial^2 F}{\partial\mu^2} = f\\\ell_y,
\qquad \frac{\partial^2 F}{\partial\mu\\\partial\sigma} =
f\left(z\ell_y + \frac{1}{\sigma}\right), \qquad \frac{\partial^2
F}{\partial\sigma^2} = f\left(z^2\ell_y + \frac{2z}{\sigma}\right).\$\$

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

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

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\\ell_y = \partial\log
f/\partial y\\.

## See also

[`loc_scale_hess_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_hess_cdf.md)
for the shared body;
[`distrib_grad_cdf.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.LogisticDistrib.md)
for the first order;
[`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

## Examples

``` r
d <- logistic_distrib()
th <- list(mu = 0.3, sigma = 1.2)
q <- c(-1, 0.5, 2)

# Against a central difference of the cdf, which shares no arithmetic.
exact <- distrib_hess_cdf(d, q, th, log = FALSE)
fd <- numerical_cdf_deriv(d, q, th, order = 2)
max(abs(unlist(exact[names(fd)]) / unlist(fd) - 1))
#> [1] 3.484743e-07
```
