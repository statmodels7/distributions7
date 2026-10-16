# Cauchy Log-CDF Gradient

Closed form: \\\partial F/\partial\mu = -f(q)\\ and \\\partial
F/\partial\sigma = -z f(q)\\ with \\z = (q-\mu)/\sigma\\, computed by
[`cauchy_cdf_method()`](https://statmodels7.github.io/distributions7/reference/cauchy_cdf_method.md)
as \\-\mathrm{Im}(\mathrm{i}^k Y^{-1})/\pi\\ with \\Y = q - \mu +
\mathrm{i}\sigma\\, without squaring the residual.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `sigma` (positive), each a
  numeric vector of length 1 or `n`.

- lower.tail:

  A single logical, `TRUE` by default; `TRUE` selects the lower tail and
  `FALSE` the upper tail.

- log:

  A single logical, `TRUE` by default; `TRUE` selects the derivatives of
  the log probability and `FALSE` those of the probability. The Cauchy's
  tails are heavy, so the probability underflows far later here than for
  a Gaussian.

## Value

A named list of two numeric vectors, `mu` and `sigma`, each the length
of `q` recycled against `theta`.

## Details

The family has no moments, and none is needed here: the distribution
function is elementary, \\F(q) = 1/2 + \arctan(z)/\pi\\, and the
identity \\\partial F/\partial\mu = -f\\ holds for every location-scale
family whether or not its moments converge. The quantities this page
returns are therefore exact where
[`mean.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.CauchyDistrib.md)
and its siblings return `NaN`.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\F\\ the distribution function.
Neither \\\mu\\ nor \\\sigma\\ is a moment: they are the median and the
half-interquartile range.

## See also

[`cauchy_cdf_method()`](https://statmodels7.github.io/distributions7/reference/cauchy_cdf_method.md)
for the body;
[`distrib_hess_cdf.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.CauchyDistrib.md)
for the second order;
[`mean.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.CauchyDistrib.md)
for the moments, which do not exist;
[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

## Examples

``` r
d <- cauchy_distrib()
th <- list(mu = 0.3, sigma = 1.2)
q <- c(-1, 0.5, 2)

# Exact, although the family has no moments.
all.equal(distrib_grad_cdf(d, q, th, log = FALSE)$mu,
          -distrib_pdf(d, q, th))
#> [1] TRUE

# The heavy tail keeps the log-scale gradient finite far out.
distrib_grad_cdf(d, -1000, th)$mu
#> [1] -0.0009996991
```
