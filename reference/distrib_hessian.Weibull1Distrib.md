# Weibull Observed Hessian

Computes the three distinct second derivatives of the Weibull
log-density with respect to the scale \\\mu\\ and the shape \\\sigma\\,
one value per observation, in closed form. With \\z = y/\mu\\ and \\u =
z^{\sigma}\\, \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
\dfrac{\sigma}{\mu^2}\left\\1 - (1 + \sigma) u\right\\, \qquad
\dfrac{\partial^2 \ell}{\partial \sigma^2} = -\dfrac{1}{\sigma^2} - u
(\log z)^2,\$\$ \$\$\dfrac{\partial^2 \ell}{\partial \mu \\ \partial
\sigma} = \dfrac{1}{\mu}\left(u - 1 + \sigma u \log z\right).\$\$ All
three carry the data through \\u\\, so none is free of \\y\\ and the
observed matrix differs from its expectation at every observation.

## Arguments

- distrib:

  A `Weibull1Distrib` object, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- y:

  A numeric vector of positive observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of three numeric vectors, `mu_mu`, `sigma_sigma` and
`mu_sigma`, each of length `max(length(y), length(mu), length(sigma))`.
The three name the distinct entries of a symmetric \\2 \times 2\\ matrix
per observation.

## See also

[`distrib_gradient.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Weibull1Distrib.md)
for the score,
[`distrib_expected_hessian.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Weibull1Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Weibull1Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- weibull1_distrib()
y <- c(0.5, 1.2, 3.0)
th <- list(mu = 2, sigma = 1.5)
h <- distrib_hessian(d, y, th)

# The three closed forms, written out.
u <- (y / 2)^1.5; lz <- log(y / 2)
all.equal(h$mu_mu, 1.5 * (1 - 2.5 * u) / 4)
#> [1] TRUE
all.equal(h$sigma_sigma, -1 / 1.5^2 - u * lz^2)
#> [1] TRUE
all.equal(h$mu_sigma, (u - 1 + 1.5 * u * lz) / 2)
#> [1] TRUE

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-5
up <- distrib_gradient(d, y, list(mu = 2 + eps, sigma = 1.5))$mu
dn <- distrib_gradient(d, y, list(mu = 2 - eps, sigma = 1.5))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
#> [1] TRUE

# The curvature in mu is positive wherever u < 1/(1 + sigma), that is
# below the 33rd percentile at this shape, so the observed information is
# not positive definite at every observation while its expectation is.
h$mu_mu
#> [1]  0.25781250 -0.06071063 -1.34729748
```
