# Logistic Log-CDF Gradient

Closed form, from the location-scale structure: \\\partial F/\partial\mu
= -f(q)\\ and \\\partial F/\partial\sigma = -z f(q)\\ with \\z =
(q-\mu)/\sigma\\. The method is
[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
itself. For this family the density is \\F(1-F)/\sigma\\, so the
log-scale gradient is elementary in \\F\\ alone and no special
evaluation is needed anywhere on the line.

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

A named list of two numeric vectors, `mu` and `sigma`, each the length
of `q` recycled against `theta`.

## Notation

\\\mu\\ is the location, \\\sigma \> 0\\ the scale, \\z =
(q-\mu)/\sigma\\, \\f\\ the density and \\F\\ the distribution function.

## See also

[`loc_scale_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/loc_scale_grad_cdf.md)
for the shared body;
[`distrib_hess_cdf.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hess_cdf.LogisticDistrib.md)
for the second order;
[`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

## Examples

``` r
d <- logistic_distrib()
th <- list(mu = 0.3, sigma = 1.2)
q <- c(-1, 0.5, 2)

# The mean component on the log scale is -f/F, which here is -(1-F)/sigma.
Fq <- distrib_cdf(d, q, th)
all.equal(distrib_grad_cdf(d, q, th)$mu, -(1 - Fq) / 1.2)
#> [1] TRUE
```
