# Weibull Third-Order Derivatives

Computes the four distinct third derivatives of the Weibull log-density
in \\\mu\\ and \\\sigma\\, in closed form, observed at the data or
expected under the model. With \\z = y/\mu\\, \\u = z^{\sigma}\\ and \\L
= \log z\\, each component is a polynomial in \\u\\, \\uL\\, \\uL^2\\
and \\uL^3\\ with coefficients rational in \\\mu\\ and \\\sigma\\; the
pure-\\\mu\\ component, for instance, is \$\$\dfrac{\partial^3
\ell}{\partial \mu^3} = \dfrac{\sigma}{\mu^3}\left\\-2 +
(1+\sigma)(2+\sigma) u\right\\.\$\$ The expected values replace each
\\\mathbb{E}\[u L^k\]\\ by \\\Gamma^{(k)}(2)/\sigma^k\\, so
`expected = TRUE` is exact here and needs no quadrature. The arithmetic
runs in a compiled kernel decomposed over the elements of the output, so
the result does not depend on the thread count.

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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
`mu_sigma_sigma` and `sigma_sigma_sigma`, each of length
`max(length(y), length(mu), length(sigma))`. The four name the distinct
entries of a symmetric third-order array over two parameters.

## See also

[`distrib_hessian.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Weibull1Distrib.md)
for the order below,
[`distrib_deriv4.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Weibull1Distrib.md)
for the order above, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- weibull1_distrib()
y <- c(0.5, 1.2, 3.0)
th <- list(mu = 2, sigma = 1.5)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_sigma_sigma"   
#> [4] "sigma_sigma_sigma"

# The pure-mu component, written out.
u <- (y / 2)^1.5
all.equal(d3$mu_mu_mu, 1.5 * (-2 + 2.5 * 3.5 * u) / 8)
#> [1] TRUE

# A central difference of the Hessian reproduces the mixed component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 2 + eps, sigma = 1.5))$mu_sigma
dn <- distrib_hessian(d, y, list(mu = 2 - eps, sigma = 1.5))$mu_sigma
all.equal((up - dn) / (2 * eps), d3$mu_mu_sigma, tolerance = 1e-6)
#> [1] TRUE

# The expected branch is exact, and averaging the observed one over draws
# reaches it.
set.seed(6)
z <- distrib_rng(d, 2e5, th)
rbind(expected = vapply(distrib_deriv3(d, y, th, expected = TRUE),
                        function(v) v[1], numeric(1)),
      averaged = vapply(distrib_deriv3(d, z, th), mean, numeric(1)))
#>          mu_mu_mu mu_mu_sigma mu_sigma_sigma sigma_sigma_sigma
#> expected 1.265625   -1.014240      0.5564164         0.4475670
#> averaged 1.265107   -1.014251      0.5574974         0.4452507
```
