# Gumbel Third-Order Derivatives

Computes the four distinct third derivatives of the Gumbel log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form. With \\z =
(y-\mu)/\sigma\\ and \\w = e^{-z}\\, every component is a polynomial in
\\z\\ and in \\z^j w\\ divided by a power of \\\sigma\\.

With `expected = TRUE` the expectations are returned, also in closed
form. They rest on \\w\\ being standard exponential under the model,
which makes \\\mathbb{E}\[z^k w\] = (-1)^k \Gamma^{(k)}(2)\\ and
\\\mathbb{E}\[z\] = \gamma\\: every expectation the family needs is a
derivative of \\\Gamma\\ at 2, assembled from polygamma functions there.
Since both routes are closed form, `approx` and `nsim` are ignored.

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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
`mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`. The names enumerate the
distinct multi-indices of order three in two parameters.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-density in
parameters \\i\\, \\j\\ and \\k\\. \\z = (y-\mu)/\sigma\\, \\w =
e^{-z}\\ and \\\gamma\\ is the Euler-Mascheroni constant.

## See also

[`distrib_hessian.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GumbelDistrib.md)
for the order below and
[`distrib_deriv4.GumbelDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GumbelDistrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- gumbel_distrib()
y <- c(-1, 0, 1)
th <- list(mu = 0, sigma = 1)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_sigma_sigma"   
#> [4] "sigma_sigma_sigma"

# The pure location component is -w/sigma^3, which is the curvature in mu
# divided by sigma. Shown at a scale where the two differ.
th2 <- list(mu = 0, sigma = 2)
all.equal(distrib_deriv3(d, y, th2)$mu_mu_mu,
          distrib_hessian(d, y, th2)$mu_mu / 2)
#> [1] TRUE

# The expected values are constants at a fixed parameter setting.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] -1
#> 
#> $mu_mu_sigma
#> [1] 2.422784
#> 
#> $mu_sigma_sigma
#> [1] -2.514818
#> 
#> $sigma_sigma_sigma
#> [1] 9.431545
#> 

# A central difference of the Hessian reproduces the same component.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 0 + eps, sigma = 1))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 0 - eps, sigma = 1))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-5)
#> [1] TRUE
```
