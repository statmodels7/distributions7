# Inverse Gaussian Observed Hessian in Mean and Dispersion

Computes the three distinct second derivatives of the inverse Gaussian
log-density with respect to \\\mu\\ and \\\phi\\, one value per
observation, in closed form: \$\$\ell^{(\mu\mu)} = -\dfrac{3y -
2\mu}{\phi\mu^4}, \qquad \ell^{(\phi\phi)} = \dfrac{\phi -
2(y-\mu)^2/(\mu^2 y)}{2\phi^3}, \qquad \ell^{(\mu\phi)} = -\dfrac{y -
\mu}{\phi^2\mu^3}.\$\$

**Neither diagonal entry is negative at every observation.** The
curvature in \\\mu\\ turns positive wherever \\y \< 2\mu/3\\, and the
curvature in \\\phi\\ is negative only where \\\phi \< 2(y-\mu)^2/(\mu^2
y)\\, so at \\y = \mu\\ it is \\+1/(2\phi^2)\\. An observation near the
mean therefore contributes positive curvature in the dispersion, and a
Newton step taken on the observed Hessian can move the wrong way. The
expected Hessian is negative definite everywhere, which is why
[`fisher_scoring()`](https://statmodels7.github.io/distributions7/reference/fisher_scoring.md)
is the more stable route here; see
[`distrib_expected_hessian.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.InvGauss1Distrib.md).

## Arguments

- distrib:

  An `InvGauss1Distrib` object, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- y:

  A numeric vector of strictly positive observations.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  Both must be strictly positive.

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

A named list of three numeric vectors, `mu_mu`, `phi_phi` and `mu_phi`,
in that order, each of length `max(length(y), length(mu), length(phi))`.
The three name the distinct entries of a symmetric \\2 \times 2\\ matrix
per observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-density in
parameters \\i\\ and \\j\\; parenthesized superscripts name derivatives.

## See also

[`distrib_gradient.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.InvGauss1Distrib.md)
for the score,
[`distrib_expected_hessian.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.InvGauss1Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.InvGauss1Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- invgauss1_distrib()
y <- c(0.5, 1, 2)
th <- list(mu = 1, phi = 2)
h <- distrib_hessian(d, y, th)
h
#> $mu_mu
#> [1]  0.25 -0.50 -2.00
#> 
#> $phi_phi
#> [1] 0.0625 0.1250 0.0625
#> 
#> $mu_phi
#> [1]  0.125  0.000 -0.250
#> 

# The three closed forms, written out.
all.equal(h$mu_mu, -(3 * y - 2 * 1) / (2 * 1^4))
#> [1] TRUE
all.equal(h$phi_phi, (2 - 2 * (y - 1)^2 / (1^2 * y)) / (2 * 2^3))
#> [1] TRUE
all.equal(h$mu_phi, -(y - 1) / (2^2 * 1^3))
#> [1] TRUE

# Positive curvature in mu below 2 mu/3, and in phi at y = mu.
c(mu_mu_at_half = distrib_hessian(d, 0.5, th)$mu_mu,
  phi_phi_at_mu = distrib_hessian(d, 1, th)$phi_phi,
  one_over_2phi2 = 1 / (2 * 2^2))
#>  mu_mu_at_half  phi_phi_at_mu one_over_2phi2 
#>          0.250          0.125          0.125 

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 1, phi = 2 + eps))$mu
dn <- distrib_gradient(d, y, list(mu = 1, phi = 2 - eps))$mu
all.equal((up - dn) / (2 * eps), h$mu_phi, tolerance = 1e-5)
#> [1] TRUE
```
