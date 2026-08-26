# Exponential Third-Order Derivative

Computes the third derivative of the exponential log-density with
respect to the mean, in closed form: \$\$\dfrac{\partial^3
\ell}{\partial \mu^3} = -\dfrac{2}{\mu^3} + \dfrac{6y}{\mu^4}.\$\$ With
`expected = TRUE` the expectation is returned, obtained by substituting
\\\mathbb{E}\[Y\] = \mu\\, which gives \\4/\mu^3\\. Both routes are
closed form, so no quadrature is run and `approx` and `nsim` are
ignored.

The family has one parameter, so there is one component, where a
two-parameter family carries four at this order.

## Arguments

- distrib:

  An `ExponentialDistrib` object, from
  [`exponential_distrib()`](https://statmodels7.github.io/distributions7/reference/exponential_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length is used.

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

A named list of one numeric vector, `mu_mu_mu`, of length
`max(length(y), length(mu))`.

## Notation

\\\ell^{(\mu\mu\mu)}\\ is the third derivative of the log-density with
respect to \\\mu\\. Parenthesized superscripts name derivatives; a
subscript on \\\ell\\ never does.

## See also

[`distrib_hessian.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.ExponentialDistrib.md)
for the order below and
[`distrib_deriv4.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.ExponentialDistrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- exponential_distrib()
y <- c(0.3, 1.1, 4.0)
th <- list(mu = 2)

# The closed form, written out.
all.equal(distrib_deriv3(d, y, th)$mu_mu_mu, -2 / 2^3 + 6 * y / 2^4)
#> [1] TRUE

# The expected value is 4/mu^3.
unique(distrib_deriv3(d, y, th, expected = TRUE)$mu_mu_mu)
#> [1] 0.5
4 / 2^3
#> [1] 0.5

# A central difference of the Hessian reproduces the observed value.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 2 + eps))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 2 - eps))$mu_mu
all.equal((up - dn) / (2 * eps), distrib_deriv3(d, y, th)$mu_mu_mu,
          tolerance = 1e-5)
#> [1] TRUE
```
