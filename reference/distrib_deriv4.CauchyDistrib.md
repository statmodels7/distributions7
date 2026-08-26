# Cauchy Fourth-Order Derivatives

Computes the five distinct fourth derivatives of the Cauchy log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form. As at third
order, every component is a rational function of \\r = y - \mu\\ and
\\\sigma\\ over a power of \\d = \sigma^2 + r^2\\ and so is bounded in
\\y\\.

With `expected = TRUE` the expectations under the model are returned,
also in closed form; the two components odd in \\r\\ vanish by symmetry.
Both routes are closed form, so `approx` and `nsim` are ignored.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

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

\\\ell^{(i j k l)}\\ is the fourth derivative of the log-density with
respect to parameters \\i\\, \\j\\, \\k\\ and \\l\\. Parenthesized
superscripts name derivatives.

## See also

[`distrib_deriv3.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.CauchyDistrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- cauchy_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
names(distrib_deriv4(d, y, th))
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_sigma"         
#> [3] "mu_mu_sigma_sigma"       "mu_sigma_sigma_sigma"   
#> [5] "sigma_sigma_sigma_sigma"

# Expected values: the two components odd in the residual vanish.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] 0.1481481
#> 
#> $mu_mu_mu_sigma
#> [1] 0
#> 
#> $mu_mu_sigma_sigma
#> [1] -0.1481481
#> 
#> $mu_sigma_sigma_sigma
#> [1] 0
#> 
#> $sigma_sigma_sigma_sigma
#> [1] -1.037037
#> 

# A central difference of the third order reproduces the fourth.
eps <- 1e-5
up <- distrib_deriv3(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu_mu_mu
all.equal((up - dn) / (2 * eps), distrib_deriv4(d, y, th)$mu_mu_mu_mu,
          tolerance = 1e-5)
#> [1] TRUE
```
