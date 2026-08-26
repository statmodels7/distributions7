# Poisson Fourth-Order Derivative

Computes the fourth derivative of the Poisson log-mass with respect to
the mean, in closed form: \$\$\dfrac{\partial^4 \ell}{\partial \mu^4} =
-\dfrac{6y}{\mu^4}.\$\$ With `expected = TRUE` the expectation is
returned, obtained by substituting \\\mathbb{E}\[Y\] = \mu\\, which
gives \\-6/\mu^3\\. Both routes are closed form, so `approx` and `nsim`
are ignored.

## Arguments

- distrib:

  A `PoissonDistrib` object, from
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

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

  Ignored, a closed form being available for both the observed and the
  expected value.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of one numeric vector, `mu_mu_mu_mu`, of length
`max(length(y), length(mu))`.

## Notation

\\\ell^{(\mu\mu\mu\mu)}\\ is the fourth derivative of the log-mass with
respect to \\\mu\\. Parenthesized superscripts name derivatives.

## See also

[`distrib_deriv3.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.PoissonDistrib.md)
for the order below;
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- poisson_distrib()
y <- c(0, 2, 7)
th <- list(mu = 3)

# The closed form, written out.
all.equal(distrib_deriv4(d, y, th)$mu_mu_mu_mu, -6 * y / 3^4)
#> [1] TRUE

# The expected value is -6/mu^3.
unique(distrib_deriv4(d, y, th, expected = TRUE)$mu_mu_mu_mu)
#> [1] -0.2222222
-6 / 3^3
#> [1] -0.2222222

# A central difference of the third order reproduces the fourth.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(mu = 3 + eps))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 3 - eps))$mu_mu_mu
all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
          tolerance = 1e-5)
#> [1] TRUE
```
