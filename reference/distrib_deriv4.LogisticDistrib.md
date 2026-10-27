# Logistic Fourth-Order Derivatives

Computes the five distinct fourth derivatives of the logistic
log-density with respect to \\\mu\\ and \\\sigma\\. The observed values
are closed form, as polynomials in \\z = (y-\mu)/\sigma\\ with the
\\g_j\\ of
[`distrib_deriv3.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LogisticDistrib.md)
as coefficients: \$\$\dfrac{\partial^4 \ell}{\partial \mu^4} =
\dfrac{g_4}{\sigma^4}, \qquad \dfrac{\partial^4 \ell}{\partial \mu^3
\partial \sigma} = \dfrac{3g_3 + z g_4}{\sigma^4}, \qquad
\dfrac{\partial^4 \ell}{\partial \mu^2 \partial \sigma^2} =
\dfrac{6g_2 + 6z g_3 + z^2 g_4}{\sigma^4},\$\$ \$\$\dfrac{\partial^4
\ell}{\partial \mu \partial \sigma^3} = \dfrac{6g_1 + 18z g_2 + 9z^2
g_3 + z^3 g_4}{\sigma^4}, \qquad \dfrac{\partial^4 \ell}{\partial
\sigma^4} = \dfrac{6 + 24z g_1 + 36z^2 g_2 + 12z^3 g_3 + z^4
g_4}{\sigma^4}.\$\$

With `expected = TRUE` the expectations under the model are returned in
closed form, as at third order: \\1/(15\sigma^4)\\ for `mu_mu_mu_mu`,
\\(\pi^2-60)/(45\sigma^4)\\ for `mu_mu_sigma_sigma`,
\\(7\pi^4-600\pi^2-810)/(225\sigma^4)\\ for `sigma_sigma_sigma_sigma`,
and zero for the two components odd in \\z\\.

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectations under the model are
  returned, in closed form. Defaults to `FALSE`.

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

  A single positive integer, how many threads the kernel may use for the
  observed values. Defaults to `1L`.

## Value

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_sigma`,
`mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
`sigma_sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell^{(i j k l)}\\ is the fourth derivative of the log-density with
respect to parameters \\i\\, \\j\\, \\k\\ and \\l\\. Parenthesized
superscripts name derivatives.

## See also

[`distrib_deriv3.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LogisticDistrib.md)
for the order below and for the \\g_j\\ used here;
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
names(distrib_deriv4(d, y, th))
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_sigma"         
#> [3] "mu_mu_sigma_sigma"       "mu_sigma_sigma_sigma"   
#> [5] "sigma_sigma_sigma_sigma"

# A central difference of the third order reproduces the fourth.
eps <- 1e-5
up <- distrib_deriv3(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu_mu_mu
all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
          tolerance = 1e-5)
#> [1] TRUE
```
