# Gumbel Log-CDF Gradient

Closed form from the location-scale structure, as for the Gaussian:
\\\partial F/\partial\mu = -f(q)\\ and \\\partial F/\partial\sigma = -z
f(q)\\ with \\z = (q-\mu)/\sigma\\. The method is
[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
itself.

## Arguments

- distrib:

  A `GumbelDistrib` object, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

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

A named list of two numeric vectors, `mu` and `sigma`, each the length
of `q` recycled against `theta`.

## Details

The family's distribution function is \\F(q) = \exp(-e^{-z})\\, so
nothing here needs a quadrature or a series, which is why a Gumbel is
usable for censored extreme-value data: the score of a right-censored
observation is \\f(q)/S(q)\\ and comes from the density alone.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\F\\ the distribution function.
The mean is \\\mu + \gamma\sigma\\ and not \\\mu\\.

## See also

[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
for the shared body;
[`distrib_hess_cdf.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.GumbelDistrib.md)
for the second order;
[`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

## Examples

``` r
d <- gumbel_distrib()
th <- list(mu = 0.3, sigma = 1.2)
q <- c(-1, 0.5, 2)

all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
          -distrib_pdf(d, q, th))
#> [1] TRUE
```
