# Gumbel Fourth-Order Derivatives

Computes the five distinct fourth derivatives of the Gumbel log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form, in the notation
of
[`distrib_deriv3.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GumbelDistrib.md):
with \\z = (y-\mu)/\sigma\\ and \\w = e^{-z}\\, every component is a
polynomial in \\z\\ and in \\z^j w\\ divided by a power of \\\sigma\\.

With `expected = TRUE` the expectations are returned, also closed form,
resting on the same fact: \\w\\ is standard exponential under the model,
so \\\mathbb{E}\[z^k w\] = (-1)^k \Gamma^{(k)}(2)\\. `approx` and `nsim`
are ignored either way.

## Arguments

- distrib:

  A `GumbelDistrib` object, from
  [`gumbel_distrib()`](https://statmodels7.github.io/distributions7/reference/gumbel_distrib.md).

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

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_sigma`,
`mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
`sigma_sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density in
parameters \\i\\, \\j\\, \\k\\ and \\l\\. \\z = (y-\mu)/\sigma\\ and \\w
= e^{-z}\\.

## See also

[`distrib_deriv3.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.GumbelDistrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic, and
[`distrib_expected_hessian.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.GumbelDistrib.md)
for the second-order expectation these extend.

## Examples

``` r
d <- gumbel_distrib()
y <- c(-1, 0, 1)
th <- list(mu = 0, sigma = 1)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_sigma"         
#> [3] "mu_mu_sigma_sigma"       "mu_sigma_sigma_sigma"   
#> [5] "sigma_sigma_sigma_sigma"

# The expected values are constants at a fixed parameter setting.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] -1
#> 
#> $mu_mu_mu_sigma
#> [1] 3.422784
#> 
#> $mu_mu_sigma_sigma
#> [1] -9.360387
#> 
#> $mu_sigma_sigma_sigma
#> [1] 15.51271
#> 
#> $sigma_sigma_sigma_sigma
#> [1] -55.30802
#> 

# A central difference of the third order reproduces it.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(mu = 0 + eps, sigma = 1))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 0 - eps, sigma = 1))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-4)
#> [1] TRUE
```
