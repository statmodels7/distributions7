# Geometric Third-Order Derivative

Computes the third derivative of the geometric log-mass with respect to
the mean, in closed form. The log-mass is \\y\log\mu -
(y+1)\log(1+\mu)\\, so every derivative is a combination of \\\mu^{-k}\\
and \\(1+\mu)^{-k}\\, and the data enter linearly through \\y\\. The
value is \\2(y/\mu^3 - (y+1)/(1+\mu)^3)\\.

The kernel writes the bracket \\y/\mu^k - (y+1)/(1+\mu)^k\\ of order
\\k\\ as \\y(1 - q^k)/\mu^k - 1/(1+\mu)^k\\ with \\q = \mu/(1+\mu)\\,
and computes \\1 - q^k\\ as
\\-\mathrm{expm1}(-k\\\mathrm{log1p}(1/\mu))\\. The two terms of the
original form agree to a relative \\1/\mu\\ where \\y\\ is near \\\mu\\,
and their difference would lose that many digits at a large mean; in
this form the value keeps full relative accuracy. A term with a zero
coefficient contributes zero, so the value at \\y = 0\\ stays finite
however small \\\mu\\ is.

With `expected = TRUE` the expectation is returned, obtained by
substituting \\\mathbb{E}\[Y\] = \mu\\, which gives \\2(1 - q^2)/\mu^2\\
with the same device. Both routes are closed form and the method runs no
quadrature.

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
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic reads it, and this method does not.

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  The generic applies an object before dispatching on the family, so
  this method receives `NULL`, and it returns the closed form of the
  family. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

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

# The expected value is one number because the data term is averaged out.
unique(distrib_deriv3(d, y, th, expected = TRUE)$mu_mu_mu)
#> [1] 0.09722222
```
