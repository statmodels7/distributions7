# Cauchy Third-Order Derivatives

Computes the four distinct third derivatives of the Cauchy log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form. Every component
is a rational function of \\r = y - \mu\\ and \\\sigma\\ over a power of
\\d = \sigma^2 + r^2\\, so all of them are bounded in \\y\\ and decay as
\\\|r\|\\ grows.

With `expected = TRUE` the expectations under the model are returned,
also in closed form: the two components odd in \\r\\ vanish by symmetry,
and the other two are \\1/(2\sigma^3)\\ and \\3/(2\sigma^3)\\. Both
routes are closed form, so no quadrature is run and `approx` and `nsim`
are ignored.

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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
`mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`.

## Notation

\\\ell^{(i j k)}\\ is the third derivative of the log-density with
respect to parameters \\i\\, \\j\\ and \\k\\. Parenthesized superscripts
name derivatives; a subscript on \\\ell\\ never does.

## See also

[`distrib_hessian.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.CauchyDistrib.md)
for the order below and
[`distrib_deriv4.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.CauchyDistrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- cauchy_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_sigma_sigma"   
#> [4] "sigma_sigma_sigma"

# Expected values: the two components odd in the residual vanish.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] 0
#> 
#> $mu_mu_sigma
#> [1] 0.1481481
#> 
#> $mu_sigma_sigma
#> [1] 0
#> 
#> $sigma_sigma_sigma
#> [1] 0.4444444
#> 
c(1 / (2 * 1.5^3), 3 / (2 * 1.5^3))
#> [1] 0.1481481 0.4444444

# A central difference of the Hessian reproduces the observed component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#> [1] TRUE
```
