# Weibull Fourth-Order Derivatives

Computes the five distinct fourth derivatives of the Weibull log-density
in \\\mu\\ and \\\sigma\\, in closed form, observed at the data or
expected under the model. The construction is the one
[`distrib_deriv3.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Weibull1Distrib.md)
describes carried one order further: with \\z = y/\mu\\, \\u =
z^{\sigma}\\ and \\L = \log z\\, each component is a polynomial in \\u
L^k\\ for \\k \le 4\\, and the pure-\\\mu\\ component is
\$\$\dfrac{\partial^4 \ell}{\partial \mu^4} =
\dfrac{\sigma}{\mu^4}\left\\6 - (1+\sigma)(2+\sigma)(3+\sigma)
u\right\\.\$\$ The expected values replace each \\\mathbb{E}\[u L^k\]\\
by \\\Gamma^{(k)}(2)/\sigma^k\\, so `expected = TRUE` is exact and needs
no quadrature. The arithmetic runs in a compiled kernel decomposed over
the elements of the output, so the result does not depend on the thread
count.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- y:

  A numeric vector of positive observations. With `expected = TRUE` only
  its length is read.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectation under the model is
  returned in place of the value at the data. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored here, both branches being exact. Accepted so that the
  signature matches the generic's.

- nsim:

  Ignored here, for the same reason. Defaults to `10000`.

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

## See also

[`distrib_deriv3.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Weibull1Distrib.md)
for the order below,
[`distrib_hessian.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Weibull1Distrib.md)
for the second order, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- weibull1_distrib()
y <- c(0.5, 1.2, 3.0)
th <- list(mu = 2, sigma = 1.5)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_sigma"         
#> [3] "mu_mu_sigma_sigma"       "mu_sigma_sigma_sigma"   
#> [5] "sigma_sigma_sigma_sigma"

# The pure-mu component, written out.
u <- (y / 2)^1.5
all.equal(d4$mu_mu_mu_mu, 1.5 * (6 - 2.5 * 3.5 * 4.5 * u) / 16)
#> [1] TRUE

# A central difference of the third order reproduces it.
eps <- 1e-4
up <- distrib_deriv3(d, y, list(mu = 2 + eps, sigma = 1.5))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 2 - eps, sigma = 1.5))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-5)
#> [1] TRUE

# The expected branch is exact, and averaging the observed one over draws
# reaches it.
set.seed(8)
z <- distrib_rng(d, 2e5, th)
rbind(expected = vapply(distrib_deriv4(d, y, th, expected = TRUE),
                        function(v) v[1], numeric(1)),
      averaged = vapply(distrib_deriv4(d, z, th), mean, numeric(1)))
#>          mu_mu_mu_mu mu_mu_mu_sigma mu_mu_sigma_sigma mu_sigma_sigma_sigma
#> expected   -3.128906       2.431170         -1.406913            0.6578897
#> averaged   -3.131105       2.434947         -1.412335            0.6632885
#>          sigma_sigma_sigma_sigma
#> expected               -1.537180
#> averaged               -1.540637
```
