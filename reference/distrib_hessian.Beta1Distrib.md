# Beta Observed Hessian in Mean and Precision

Computes the three distinct second derivatives of the beta log-density
with respect to \\\mu\\ and \\\phi\\, one value per observation, in
closed form. With \\\alpha = \mu\phi\\, \\\beta = (1-\mu)\phi\\ and
\\\psi_1\\ the trigamma function, \$\$\ell^{(\mu\mu)} =
-\phi^2\left\\\psi_1(\alpha) + \psi_1(\beta)\right\\, \qquad
\ell^{(\phi\phi)} = \psi_1(\phi) - \mu^2\psi_1(\alpha) -
(1-\mu)^2\psi_1(\beta),\$\$ \$\$\ell^{(\mu\phi)} = \log\dfrac{y}{1-y} -
\psi(\alpha) + \psi(\beta) - \phi\left\\\mu\psi_1(\alpha) -
(1-\mu)\psi_1(\beta)\right\\.\$\$

Two of the three are free of the data. The mixed entry is not: it
carries the log-odds residual \\\log\\y/(1-y)\\ - \psi(\alpha) +
\psi(\beta)\\, which is the mean component of the score divided by
\\\phi\\ and has expectation zero. So this Hessian and
[`distrib_expected_hessian.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Beta1Distrib.md)
agree in `mu_mu` and `phi_phi` and differ in `mu_phi` by exactly that
residual.

## Arguments

- distrib:

  A `Beta1Distrib` object, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- y:

  A numeric vector of observations in \\(0, 1)\\. Only the mixed entry
  reads it.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `mu` must lie strictly in \\(0, 1)\\ and `phi` must be strictly
  positive.

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
\\\psi\\ and \\\psi_1\\ are the digamma and trigamma functions.

## See also

[`distrib_gradient.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Beta1Distrib.md)
for the score,
[`distrib_expected_hessian.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Beta1Distrib.md)
for the expectation of this quantity,
[`distrib_deriv3.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Beta1Distrib.md)
for the order above, and
[`distrib_hessian()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.md)
for the generic.

## Examples

``` r
d <- beta1_distrib()
y <- c(0.2, 0.5, 0.8)
th <- list(mu = 0.4, phi = 5)
h <- distrib_hessian(d, y, th)

# Two entries are constant across the observations and one is not.
h
#> $mu_mu
#> [1] -25.9967 -25.9967 -25.9967
#> 
#> $phi_phi
#> [1] -0.02404276 -0.02404276 -0.02404276
#> 
#> $mu_phi
#> [1] -0.9913603  0.3949341  1.7812284
#> 

# The mixed entry differs from its expectation by the log-odds residual.
a <- 0.4 * 5
b <- 0.6 * 5
all.equal(h$mu_phi - distrib_expected_hessian(d, y, th)$mu_phi,
          log(y / (1 - y)) - digamma(a) + digamma(b))
#> [1] TRUE

# It is the second derivative of the log-density, so a central difference
# of the score reproduces it.
eps <- 1e-6
up <- distrib_gradient(d, y, list(mu = 0.4, phi = 5 + eps))$mu
dn <- distrib_gradient(d, y, list(mu = 0.4, phi = 5 - eps))$mu
all.equal((up - dn) / (2 * eps), h$mu_phi, tolerance = 1e-5)
#> [1] TRUE
```
