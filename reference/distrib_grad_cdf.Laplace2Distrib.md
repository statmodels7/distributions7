# Laplace Log-CDF Gradient in Location and Rate

Closed form: \\\partial F/\partial\mu = -f(q)\\ and \\\partial
F/\partial\lambda = (q-\mu)\\f(q)/\lambda\\. This is the same law as
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
written by its rate \\\lambda = 1/\sigma\\, and the chain rule turns the
scale component into the one above.

## Arguments

- distrib:

  A `Laplace2Distrib` object, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- q:

  A numeric vector of quantiles.

- theta:

  A named list with components `mu` and `lambda` (positive), each a
  numeric vector of length 1 or `n`.

- lower.tail:

  Is the lower tail wanted? A single logical, `TRUE` by default.

- log:

  Are derivatives of the log probability wanted? A single logical,
  `TRUE` by default.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of two numeric vectors, `mu` and `lambda`, each the length
of `q` recycled against `theta`.

## Notation

\\\mu\\ is the location, \\\lambda \> 0\\ the rate, \\f\\ the density
and \\F\\ the distribution function. The variance is \\2/\lambda^2\\.

## See also

[`distrib_hess_cdf.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.Laplace2Distrib.md),
where the kink shows;
[`distrib_grad_cdf.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.LaplaceDistrib.md)
for the scale parametrization;
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

## Examples

``` r
d <- laplace2_distrib()
th <- list(mu = 0.3, lambda = 1 / 1.2)
q <- c(-1, 0.5, 2)

g <- distrib_grad_cdf(d, q, th, log = FALSE)
all.equal(g$mu, -distrib_pdf(d, q, th))
#> [1] TRUE
all.equal(g$lambda, (q - 0.3) * distrib_pdf(d, q, th) / (1 / 1.2))
#> [1] TRUE
```
