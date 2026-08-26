# Inverse Gaussian Fourth-Order Derivatives in Mean and Dispersion

Computes the five distinct fourth derivatives of the inverse Gaussian
log-density with respect to \\\mu\\ and \\\phi\\, in closed form, by the
same route the lower orders take: the log-density is a linear
combination of \\y\\, \\1/y\\ and \\\log y\\ with coefficients rational
in the parameters.

With `expected = TRUE` the expectations are returned, obtained by
replacing \\Y\\ with \\\mu\\ and \\1/Y\\ with \\1/\mu + \phi\\. Both
routes are closed form, so `approx` and `nsim` are ignored.

## Arguments

- distrib:

  An `InvGauss1Distrib` object, from
  [`invgauss1_distrib()`](https://statmodels7.github.io/distributions7/reference/invgauss1_distrib.md).

- y:

  A numeric vector of strictly positive observations. With
  `expected = TRUE` only its length is used.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  Both must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectations under the model are
  returned in place of the observed values. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available for both the observed and the
  expected values.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_phi`,
`mu_mu_phi_phi`, `mu_phi_phi_phi` and `phi_phi_phi_phi`, each of length
`max(length(y), length(mu), length(phi))`.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density in
parameters \\i\\, \\j\\, \\k\\ and \\l\\. Parenthesized superscripts
name derivatives.

## See also

[`distrib_deriv3.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.InvGauss1Distrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic, and
[`distrib_expected_hessian.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.InvGauss1Distrib.md)
for the second-order expectation these extend.

## Examples

``` r
d <- invgauss1_distrib()
y <- c(0.5, 1, 2)
th <- list(mu = 1, phi = 2)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"     "mu_mu_mu_phi"    "mu_mu_phi_phi"   "mu_phi_phi_phi" 
#> [5] "phi_phi_phi_phi"

# The expected values are constants at a fixed parameter setting.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] -18
#> 
#> $mu_mu_mu_phi
#> [1] -1.5
#> 
#> $mu_mu_phi_phi
#> [1] -0.25
#> 
#> $mu_phi_phi_phi
#> [1] 0
#> 
#> $phi_phi_phi_phi
#> [1] -0.5625
#> 

# A central difference of the third order reproduces it.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(mu = 1 + eps, phi = 2))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 1 - eps, phi = 2))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-3)
#> [1] TRUE
```
