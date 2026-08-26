# Bernoulli Third-Order Derivative

Computes the third derivative of the Bernoulli log-mass with respect to
the probability, in closed form: \$\$\dfrac{\partial^3 \ell}{\partial
\mu^3} = \dfrac{2y}{\mu^3} - \dfrac{2(1-y)}{(1-\mu)^3}.\$\$ With
`expected = TRUE` the expectation is returned, obtained by substituting
\\\mathbb{E}\[Y\] = \mu\\, which gives \\2/\mu^2 - 2/(1-\mu)^2\\. Both
routes are closed form, so no quadrature is run and `approx` and `nsim`
are ignored.

The family has one parameter, so there is one component, where a
two-parameter family carries four at this order.

## Arguments

- distrib:

  A `BernoulliDistrib` object, from
  [`bernoulli_distrib()`](https://statmodels7.github.io/distributions7/reference/bernoulli_distrib.md).

- y:

  A numeric vector of zeros and ones. With `expected = TRUE` only its
  length is used.

- theta:

  A named list with the single component `mu`, a numeric vector of
  length 1 or of the length of `y`. `mu` must lie in \\(0, 1)\\.

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

\\\ell^{(\mu\mu\mu)}\\ is the third derivative of the log-mass with
respect to \\\mu\\. Parenthesized superscripts name derivatives.

## See also

[`distrib_hessian.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.BernoulliDistrib.md)
for the order below and
[`distrib_deriv4.BernoulliDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.BernoulliDistrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- bernoulli_distrib()
y <- c(0, 1, 1)
th <- list(mu = 0.3)

# The closed form, written out.
all.equal(distrib_deriv3(d, y, th)$mu_mu_mu,
          2 * y / 0.3^3 - 2 * (1 - y) / 0.7^3)
#> [1] TRUE

# The expected value, and zero at mu = 1/2 by symmetry.
unique(distrib_deriv3(d, y, th, expected = TRUE)$mu_mu_mu)
#> [1] 18.14059
2 / 0.3^2 - 2 / 0.7^2
#> [1] 18.14059
distrib_deriv3(d, 0, list(mu = 0.5), expected = TRUE)$mu_mu_mu
#> [1] 0

# A central difference of the Hessian reproduces the observed value.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 0.3 + eps))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 0.3 - eps))$mu_mu
all.equal((up - dn) / (2 * eps), distrib_deriv3(d, y, th)$mu_mu_mu,
          tolerance = 1e-5)
#> [1] TRUE
```
