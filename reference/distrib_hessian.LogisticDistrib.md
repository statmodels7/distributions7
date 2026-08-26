# Logistic Observed Hessian

Computes the three distinct second derivatives of the logistic
log-density with respect to \\\mu\\ and \\\sigma\\, one value per
observation, in closed form. Writing \\z = (y - \mu)/\sigma\\,
\$\$\dfrac{\partial^2 \ell}{\partial \mu^2} = -\dfrac{1}{2\sigma^2}
\mathrm{sech}^2\left(\dfrac{z}{2}\right),\$\$ \$\$\dfrac{\partial^2
\ell}{\partial \sigma^2} = \dfrac{1}{\sigma^2}\left\[1 -
\dfrac{z^2}{2}\\\mathrm{sech}^2\left(\dfrac{z}{2}\right) - 2 z
\tanh\left(\dfrac{z}{2}\right)\right\],\$\$ \$\$\dfrac{\partial^2
\ell}{\partial \mu \\ \partial \sigma} =
-\dfrac{1}{\sigma^2}\left\[\tanh\left(\dfrac{z}{2}\right) +
\dfrac{z}{2}\\\mathrm{sech}^2\left(\dfrac{z}{2}\right)\right\].\$\$

The curvature in \\\mu\\ is negative everywhere and vanishes as
\\\|z\|\\ grows, so a distant observation carries almost no information
about the location. The expectations are
[`distrib_expected_hessian.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LogisticDistrib.md).

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

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
`mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`.
The three name the distinct entries of a symmetric \\2 \times 2\\ matrix
per observation.

## See also

[`distrib_gradient.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LogisticDistrib.md)
for the score,
[`distrib_expected_hessian.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.LogisticDistrib.md)
for the expectation of this quantity, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
h <- distrib_hessian(d, y, th)

# The closed form for the location, written out.
z <- (y - 0.4) / 1.5
all.equal(h$mu_mu, -(1 - tanh(z / 2)^2) / (2 * 1.5^2))
#> [1] TRUE

# Concave in mu everywhere, and flattening as the residual grows.
round(distrib_hessian(d, 0.4 + c(0, 1.5, 3, 6, 12), th)$mu_mu, 6)
#> [1] -0.222222 -0.174766 -0.093328 -0.015700 -0.000298

# A central difference of the score reproduces it.
eps <- 1e-5
up <- distrib_gradient(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu
dn <- distrib_gradient(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
#> [1] TRUE
```
