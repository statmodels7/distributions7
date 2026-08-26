# Gaussian Observed Hessian

Computes the three distinct second derivatives of the Gaussian
log-density with respect to \\\mu\\ and \\\sigma\\, one value per
observation, in closed form: \$\$\dfrac{\partial^2 \ell}{\partial \mu^2}
= -\dfrac{1}{\sigma^2}, \qquad \dfrac{\partial^2 \ell}{\partial
\sigma^2} = \dfrac{\sigma^2 - 3(y - \mu)^2}{\sigma^4}, \qquad
\dfrac{\partial^2 \ell}{\partial \mu \\ \partial \sigma} = -\dfrac{2(y -
\mu)}{\sigma^3}.\$\$ Only the curvature in \\\mu\\ is free of the data;
the other two components vary with the residual, and their expectations
are
[`distrib_expected_hessian.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian1Distrib.md).

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

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

[`distrib_gradient.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian1Distrib.md)
for the score,
[`distrib_expected_hessian.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gaussian1Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gaussian1Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
h <- distrib_hessian(d, y, th)

# The curvature in mu is constant at -1/sigma^2; the other two are not.
h$mu_mu
#> [1] -0.4444444 -0.4444444 -0.4444444
all.equal(h$sigma_sigma, (1.5^2 - 3 * (y - 0.4)^2) / 1.5^4)
#> [1] TRUE
all.equal(h$mu_sigma, -2 * (y - 0.4) / 1.5^3)
#> [1] TRUE

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-5
up <- distrib_gradient(d, y, list(mu = 0.4 + eps, sigma = 1.5))$mu
dn <- distrib_gradient(d, y, list(mu = 0.4 - eps, sigma = 1.5))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
#> [1] TRUE
```
