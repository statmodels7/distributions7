# Gamma Observed Hessian in Mean and Dispersion

Computes the three distinct second derivatives of the gamma log-density
with respect to \\\mu\\ and \\\phi\\, one value per observation, in
closed form. With \\s = 1/\phi\\ and \\z = y/\mu\\, \$\$\ell^{(\mu\mu)}
= \dfrac{s(1 - 2z)}{\mu^2}, \qquad \ell^{(\mu\phi)} = \dfrac{s^2(1 -
z)}{\mu}, \qquad \ell^{(\phi\phi)} = s^4\left\\\dfrac{1}{s} -
\psi'(s)\right\\ + 2s^3\left\\\log s + 1 - \psi(s) + \log z -
z\right\\.\$\$ Every derivative in \\\phi\\ is the corresponding
derivative in \\s\\ carried across by the one-variable chain rule, with
\\s' = -s^2\\ and \\s'' = 2s^3\\, so each polygamma function is
evaluated once.

The curvature in \\\mu\\ turns **positive** wherever \\y \< \mu/2\\, so
the observed information of a gamma is not positive definite at every
data point. Its expectation is, which is one reason a gamma model is
fitted by Fisher scoring; see
[`distrib_expected_hessian.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma1Distrib.md).

## Arguments

- distrib:

  A `Gamma1Distrib` object, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

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

A named list of three numeric vectors, `mu_mu`, `mu_phi` and `phi_phi`,
each of length `max(length(y), length(mu), length(phi))`. The three name
the distinct entries of a symmetric \\2 \times 2\\ matrix per
observation.

## Notation

\\\ell^{(ij)}\\ is the second derivative of the log-density in
parameters \\i\\ and \\j\\; parenthesized superscripts name derivatives.
\\\psi\\ and \\\psi'\\ are the digamma and trigamma functions.

## See also

[`distrib_gradient.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma1Distrib.md)
for the score,
[`distrib_expected_hessian.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma1Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gamma1Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- gamma1_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, phi = 0.5)
h <- distrib_hessian(d, y, th)

# The curvature in mu, written out, and positive at y = 1 < mu/2 = 1.5.
s <- 1 / 0.5
z <- y / 3
all.equal(h$mu_mu, s * (1 - 2 * z) / 3^2)
#> [1] TRUE
h$mu_mu
#> [1]  0.07407407 -0.22222222 -0.51851852

# The mixed entry vanishes at y = mu.
distrib_hessian(d, 3, th)$mu_phi
#> [1] 0

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 3, phi = 0.5 + eps))$phi
dn <- distrib_gradient(d, y, list(mu = 3, phi = 0.5 - eps))$phi
all.equal((up - dn) / (2 * eps), h$phi_phi, tolerance = 1e-5)
#> [1] TRUE
```
