# Laplace Log-CDF Gradient

Closed form, from the location-scale structure: \\\partial F/\partial\mu
= -f(q)\\ and \\\partial F/\partial\sigma = -z f(q)\\ with \\z =
(q-\mu)/\sigma\\. The method is
[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
itself. Both are continuous everywhere, including at \\q = \mu\\: the
density has a kink there but no jump, and the first derivatives of \\F\\
read the density itself.

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

A named list of two numeric vectors, `mu` and `sigma`, each the length
of `q` recycled against `theta`.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\F\\ the distribution function.
The variance is \\2\sigma^2\\, so the scale is not a standard deviation.

## See also

[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
for the shared body;
[`distrib_hess_cdf.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.LaplaceDistrib.md),
where the kink does show;
[`distrib_grad_cdf.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.Laplace2Distrib.md)
for the rate parametrization;
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

## Examples

``` r
d <- laplace_distrib()
th <- list(mu = 0.3, sigma = 1.2)
q <- c(-1, 0.5, 2)

all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
          -distrib_pdf(d, q, th))
#> [1] TRUE

# Continuous through the kink at q = mu.
distrib_grad_cdf(d, 0.3 + c(-1e-8, 0, 1e-8), th, log = FALSE)$mu
#> [1] -0.4166667 -0.4166667 -0.4166667
```
