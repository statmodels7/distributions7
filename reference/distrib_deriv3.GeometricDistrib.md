# Geometric Third-Order Derivative

Computes the third derivative of the geometric log-mass with respect to
the mean, in closed form. The log-mass is \\y\log\mu -
(y+1)\log(1+\mu)\\, so every derivative is a combination of \\\mu^{-k}\\
and \\(1+\mu)^{-k}\\ with the data entering linearly through \\y\\.

With `expected = TRUE` the expectation is returned, obtained by
substituting \\\mathbb{E}\[Y\] = \mu\\. Both routes are closed form, so
no quadrature is run and `approx` and `nsim` are ignored.

The family has one parameter, so there is one component, where a
two-parameter family carries four at this order.

## Arguments

- distrib:

  A `GeometricDistrib` object, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

- y:

  A numeric vector of counts. With `expected = TRUE` only its length is
  used.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectation under the model is
  returned in place of the observed value. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available for both routes.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of one numeric vector, `mu_mu_mu`, of length
`max(length(y), length(mu))`.

## Notation

\\\ell^{(\mu\mu\mu)}\\ is the third derivative of the log-mass with
respect to \\\mu\\. Parenthesized superscripts name derivatives.

## See also

[`distrib_hessian.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GeometricDistrib.md)
for the order below and
[`distrib_deriv4.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GeometricDistrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- geometric_distrib()
y <- c(0, 2, 7)
th <- list(mu = 3)

# A central difference of the Hessian reproduces the observed value.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 3 + eps))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 3 - eps))$mu_mu
all.equal((up - dn) / (2 * eps), distrib_deriv3(d, y, th)$mu_mu_mu,
          tolerance = 1e-5)
#> [1] TRUE

# The expected value is one number, the data term having been averaged out.
unique(distrib_deriv3(d, y, th, expected = TRUE)$mu_mu_mu)
#> [1] 0.09722222
```
