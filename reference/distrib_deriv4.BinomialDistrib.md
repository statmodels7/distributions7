# Binomial Fourth-Order Derivative

Computes the fourth derivative of the binomial log-mass with respect to
the probability, in closed form: \$\$\dfrac{\partial^4 \ell}{\partial
\mu^4} = -\dfrac{6y}{\mu^4} - \dfrac{6(n-y)}{(1-\mu)^4}.\$\$ With
`expected = TRUE` the expectation is returned, obtained by substituting
\\\mathbb{E}\[Y\] = n\mu\\, which gives \\-6n/\mu^3 - 6n/(1-\mu)^3\\.
Both routes are closed form, so `approx` and `nsim` are ignored.

## Arguments

- distrib:

  A `BinomialDistrib` object, from
  [`binomial_distrib()`](https://statmodels7.github.io/distributions7/reference/binomial_distrib.md).
  Its `size` property supplies the number of trials.

- y:

  A numeric vector of counts of successes. With `expected = TRUE` only
  its length is used.

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

  Ignored, a closed form being available for both routes.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of one numeric vector, `mu_mu_mu_mu`, of length
`max(length(y), length(mu), length(distrib@size))`.

## Notation

\\\ell^{(\mu\mu\mu\mu)}\\ is the fourth derivative of the log-mass with
respect to \\\mu\\. Parenthesized superscripts name derivatives.

## See also

[`distrib_deriv3.BinomialDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.BinomialDistrib.md)
for the order below;
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- binomial_distrib(size = 10)
y <- c(0, 4, 10)
th <- list(mu = 0.3)

# The closed form, written out; negative at every admissible mu.
all.equal(distrib_deriv4(d, y, th)$mu_mu_mu_mu,
          -6 * y / 0.3^4 - 6 * (10 - y) / 0.7^4)
#> [1] TRUE

# The expected value.
unique(distrib_deriv4(d, y, th, expected = TRUE)$mu_mu_mu_mu)
#> [1] -2397.149
-6 * 10 / 0.3^3 - 6 * 10 / 0.7^3
#> [1] -2397.149

# A central difference of the third order reproduces the fourth.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(mu = 0.3 + eps))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 0.3 - eps))$mu_mu_mu
all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
          tolerance = 1e-5)
#> [1] TRUE
```
