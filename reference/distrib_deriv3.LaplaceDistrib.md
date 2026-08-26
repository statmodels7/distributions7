# Laplace Third-Order Derivatives

Computes the four distinct third derivatives of the Laplace log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form. The log-density
is \\-\log\sigma - \log 2 - \|r\|/\sigma\\ with \\r = y - \mu\\, so it
is linear in \\\mu\\ on each side of the kink and every component that
differentiates twice or more in \\\mu\\ is zero. What survives comes
from \\-\log\sigma\\ and from \\\|r\|/\sigma\\.

With `expected = TRUE` the expectations under the model are returned,
also in closed form, using \\\mathbb{E}\|r\| = \sigma\\ and
\\\mathbb{E}\[\mathrm{sign}(r)\] = 0\\. Both routes are closed form, so
no quadrature is run and `approx` and `nsim` are ignored.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length is used.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectations under the model are
  returned in place of the observed values. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available for both the observed and the
  expected values.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
`mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell^{(i j k)}\\ is the third derivative of the log-density with
respect to parameters \\i\\, \\j\\ and \\k\\. Parenthesized superscripts
name derivatives; a subscript on \\\ell\\ never does. These are the
derivatives that exist away from the kink at \\y = \mu\\.

## See also

[`distrib_hessian.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.LaplaceDistrib.md)
for the order below and
[`distrib_deriv4.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.LaplaceDistrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- laplace_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_sigma_sigma"   
#> [4] "sigma_sigma_sigma"

# Linear in mu on each side of the kink, so every component that
# differentiates twice in mu is zero.
c(d3$mu_mu_mu[1], d3$mu_mu_sigma[1])
#> [1] 0 0

# The expected values, in closed form.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_sigma
#> [1] 0
#> 
#> $mu_sigma_sigma
#> [1] 0
#> 
#> $sigma_sigma_sigma
#> [1] 1.185185
#> 

# A central difference of the Hessian reproduces the observed component,
# away from the kink.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 0.4, sigma = 1.5 + eps))$sigma_sigma
dn <- distrib_hessian(d, y, list(mu = 0.4, sigma = 1.5 - eps))$sigma_sigma
all.equal((up - dn) / (2 * eps), d3$sigma_sigma_sigma, tolerance = 1e-6)
#> [1] TRUE
```
