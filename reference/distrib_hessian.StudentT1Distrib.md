# Student t Observed Hessian

Computes the six distinct second derivatives of the location-scale
Student t log-density with respect to \\\mu\\, \\\sigma\\ and \\\nu\\,
one value per observation, in closed form. With \\r = y - \mu\\ and \\D
= \nu\sigma^2 + r^2\\, \$\$\dfrac{\partial^2 \ell}{\partial \mu^2} =
\dfrac{(\nu+1)\left(r^2 - \nu\sigma^2\right)}{D^2}, \qquad
\dfrac{\partial^2 \ell}{\partial \sigma^2} =
\dfrac{\nu\left\\\nu\sigma^4 - (3\nu+1)\sigma^2 r^2 -
r^4\right\\}{\sigma^2 D^2},\$\$ \$\$\dfrac{\partial^2 \ell}{\partial
\nu^2} = \dfrac{1}{4}\left\[-\psi_1\\\left(\dfrac{\nu}{2}\right) +
\psi_1\\\left(\dfrac{\nu+1}{2}\right) + \dfrac{2\left(\nu\sigma^4 +
r^4\right)}{\nu D^2}\right\],\$\$ \$\$\dfrac{\partial^2 \ell}{\partial
\mu \\ \partial \sigma} = -\dfrac{2\nu(\nu+1)\sigma r}{D^2}, \qquad
\dfrac{\partial^2 \ell}{\partial \mu \\ \partial \nu} =
\dfrac{r\left(r^2 - \sigma^2\right)}{D^2}, \qquad \dfrac{\partial^2
\ell}{\partial \sigma \\ \partial \nu} = \dfrac{r^2\left(r^2 -
\sigma^2\right)}{\sigma D^2}.\$\$

The curvature in \\\mu\\ **turns positive** wherever \\\|r\| \>
\sigma\sqrt{\nu}\\, the same point at which the score peaks, so the
observed information is indefinite in that direction at an outlying
observation while its expectation is negative definite everywhere.

The arithmetic runs in a compiled kernel decomposed over the elements of
the output, so the result does not depend on the thread count.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

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

A named list of six numeric vectors, `mu_mu`, `sigma_sigma`, `nu_nu`,
`mu_sigma`, `mu_nu` and `sigma_nu`, each of length
`max(length(y), length(mu), length(sigma), length(nu))`. The six name
the distinct entries of a symmetric \\3 \times 3\\ matrix per
observation.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the location,
\\\sigma \> 0\\ the scale and \\\nu \> 0\\ the degrees of freedom.
\\\psi_1\\ is the trigamma function,
[`trigamma()`](https://rdrr.io/r/base/Special.html) in R.

## See also

[`distrib_gradient.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentT1Distrib.md)
for the score,
[`distrib_expected_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.StudentT1Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.StudentT1Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 5)
h <- distrib_hessian(d, y, th)
names(h)
#> [1] "mu_mu"       "sigma_sigma" "nu_nu"       "mu_sigma"    "mu_nu"      
#> [6] "sigma_nu"   

# The location and mixed location-scale components, written out.
r <- y - 0.4; D <- 5 * 1.2^2 + r^2
all.equal(h$mu_mu, 6 * (r^2 - 5 * 1.2^2) / D^2)
#> [1] TRUE
all.equal(h$mu_sigma, -2 * 5 * 6 * 1.2 * r / D^2)
#> [1] TRUE

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-5
up <- distrib_gradient(d, y, list(mu = 0.4 + eps, sigma = 1.2, nu = 5))$mu
dn <- distrib_gradient(d, y, list(mu = 0.4 - eps, sigma = 1.2, nu = 5))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-6)
#> [1] TRUE

# The curvature in mu is positive beyond |r| = sigma * sqrt(nu) = 2.68.
rr <- c(0.5, 1, 2, 4, 8, 16)
rbind(residual = rr,
      mu_mu = 6 * (rr^2 - 5 * 1.2^2) / (5 * 1.2^2 + rr^2)^2)
#>                [,1]       [,2]       [,3]      [,4]       [,5]        [,6]
#> residual  0.5000000  1.0000000  2.0000000 4.0000000 8.00000000 16.00000000
#> mu_mu    -0.7513175 -0.5532421 -0.1530612 0.0980975 0.06722636  0.02154914
```
