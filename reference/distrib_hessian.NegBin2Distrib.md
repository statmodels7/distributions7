# Negative Binomial Observed Hessian, NB2

Computes the three distinct second derivatives of the negative binomial
log-mass with respect to \\\mu\\ and \\\theta\\, one value per
observation. Writing \\s = \theta + \mu\\ and \\\psi_1\\ for the
trigamma function, \$\$\ell^{(\mu\mu)} = \dfrac{y+\theta}{s^2} -
\dfrac{y}{\mu^2}, \qquad \ell^{(\mu\theta)} = \dfrac{y-\mu}{s^2}, \qquad
\ell^{(\theta\theta)} = \psi_1(y+\theta) - \psi_1(\theta) +
\dfrac{\mu}{\theta s} + \dfrac{y-\mu}{s^2}.\$\$

The dispersion entry cancels at large \\\theta\\ for the same reason the
score does, and the kernel handles it the same way. Here the collapse is
better than a series: with \\a = \theta\\, \\b = \theta + y\\ and \\c =
\theta + \mu\\ the three leading terms combine **exactly**,
\$\$-\dfrac{y}{ab} + \dfrac{\mu}{ac} + \dfrac{y-\mu}{c^2} =
\dfrac{(y-\mu)^2}{b\\c^2},\$\$ so what is computed is that quotient plus
the remainder \\\psi_1(b) - \psi_1(a) + y/(ab)\\, and only the remainder
needs a series. To leading order the result is \\\\(y-\mu)^2 -
y\\/\theta^3\\, which is the derivative of the score's leading term and
is how the two derivations check each other.

## Arguments

- distrib:

  A `NegBin2Distrib` object, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- y:

  A numeric vector of counts.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

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

A named list of three numeric vectors, `mu_mu`, `theta_theta` and
`mu_theta`, in that order, each of length
`max(length(y), length(mu), length(theta))`. The three name the distinct
entries of a symmetric \\2 \times 2\\ matrix per observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-mass in parameters
\\i\\ and \\j\\; parenthesized superscripts name derivatives. \\\psi\\
and \\\psi_1\\ are the digamma and trigamma functions.

## See also

[`distrib_gradient.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin2Distrib.md)
for the score,
[`distrib_expected_hessian.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin2Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin2Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- negbin2_distrib()
y <- c(0, 2, 6)
th <- list(mu = 4, theta = 2)
h <- distrib_hessian(d, y, th)
h
#> $mu_mu
#> [1]  0.05555556 -0.01388889 -0.15277778
#> 
#> $theta_theta
#> [1]  0.22222222 -0.08333333 -0.12290816
#> 
#> $mu_theta
#> [1] -0.11111111 -0.05555556  0.05555556
#> 

# The three closed forms, written out at a dispersion where the direct
# expression still has its digits.
s <- 2 + 4
all.equal(h$mu_mu, (y + 2) / s^2 - y / 4^2)
#> [1] TRUE
all.equal(h$mu_theta, (y - 4) / s^2)
#> [1] TRUE
all.equal(h$theta_theta,
          trigamma(y + 2) - trigamma(2) + 4 / (2 * s) + (y - 4) / s^2)
#> [1] TRUE

# The curvature in mu is positive at y = 0, so the observed information is
# not positive definite at every count.
h$mu_mu
#> [1]  0.05555556 -0.01388889 -0.15277778

# It is the second derivative of the log-mass, so a central difference of
# the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 4 + eps, theta = 2))$mu
dn <- distrib_gradient(d, y, list(mu = 4 - eps, theta = 2))$mu
all.equal((up - dn) / (2 * eps), h$mu_mu, tolerance = 1e-5)
#> [1] TRUE
```
