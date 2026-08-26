# Laplace Fourth-Order Derivatives

Computes the five distinct fourth derivatives of the Laplace log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form. As at third
order, the log-density is linear in \\\mu\\ on each side of the kink, so
every component differentiating twice or more in \\\mu\\ is zero and
what survives comes from the \\-\log\sigma\\ term and from
\\\|r\|/\sigma\\.

With `expected = TRUE` the expectations under the model are returned,
also in closed form. Both routes are closed form, so `approx` and `nsim`
are ignored.

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

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_sigma`,
`mu_mu_sigma_sigma`, `mu_sigma_sigma_sigma` and
`sigma_sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell^{(i j k l)}\\ is the fourth derivative of the log-density with
respect to parameters \\i\\, \\j\\, \\k\\ and \\l\\. Parenthesized
superscripts name derivatives.

## See also

[`distrib_deriv3.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.LaplaceDistrib.md)
for the order below;
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- laplace_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
names(distrib_deriv4(d, y, th))
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_sigma"         
#> [3] "mu_mu_sigma_sigma"       "mu_sigma_sigma_sigma"   
#> [5] "sigma_sigma_sigma_sigma"

# The pure-sigma component is 6/sigma^4 from the -log(sigma) term, less the
# contribution of |r|/sigma.
distrib_deriv4(d, y, th)$sigma_sigma_sigma_sigma
#> [1] -3.8716049  0.8691358 -5.4518519

# The expected values, in closed form.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_mu_sigma
#> [1] 0
#> 
#> $mu_mu_sigma_sigma
#> [1] 0
#> 
#> $mu_sigma_sigma_sigma
#> [1] 0
#> 
#> $sigma_sigma_sigma_sigma
#> [1] -3.555556
#> 

# A central difference of the third order reproduces the fourth.
eps <- 1e-5
up <- distrib_deriv3(d, y, list(mu = 0.4, sigma = 1.5 + eps))$sigma_sigma_sigma
dn <- distrib_deriv3(d, y, list(mu = 0.4, sigma = 1.5 - eps))$sigma_sigma_sigma
all.equal((up - dn) / (2 * eps),
          distrib_deriv4(d, y, th)$sigma_sigma_sigma_sigma, tolerance = 1e-5)
#> [1] TRUE
```
