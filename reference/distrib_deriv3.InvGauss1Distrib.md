# Inverse Gaussian Third-Order Derivatives in Mean and Dispersion

Computes the four distinct third derivatives of the inverse Gaussian
log-density with respect to \\\mu\\ and \\\phi\\, in closed form. The
log-density is a linear combination of \\y\\, \\1/y\\ and \\\log y\\
with coefficients rational in the parameters, so every derivative is a
rational function of \\(\mu, \phi)\\ times one of those three
statistics.

With `expected = TRUE` the expectations are returned, obtained by
replacing \\Y\\ with \\\mu\\ and \\1/Y\\ with \\1/\mu + \phi\\. Both
routes are closed form, so no quadrature is run and `approx` and `nsim`
are ignored.

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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_phi`,
`mu_phi_phi` and `phi_phi_phi`, each of length
`max(length(y), length(mu), length(phi))`. The names enumerate the
distinct multi-indices of order three in two parameters.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-density in
parameters \\i\\, \\j\\ and \\k\\. Parenthesized superscripts name
derivatives.

## See also

[`distrib_hessian.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.InvGauss1Distrib.md)
for the order below and
[`distrib_deriv4.InvGauss1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.InvGauss1Distrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic and for the numerical route a family without a closed
form takes.

## Examples

``` r
d <- invgauss1_distrib()
y <- c(0.5, 1, 2)
th <- list(mu = 1, phi = 2)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"    "mu_mu_phi"   "mu_phi_phi"  "phi_phi_phi"

# The expected values are constants at a fixed parameter setting, and the
# mixed mu-phi-phi component is exactly zero, the score in mu being linear
# in y and the two parameters entering it as a product.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] 3
#> 
#> $mu_mu_phi
#> [1] 0.25
#> 
#> $mu_phi_phi
#> [1] 0
#> 
#> $phi_phi_phi
#> [1] 0.25
#> 

# A central difference of the Hessian reproduces the same component.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 1 + eps, phi = 2))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 1 - eps, phi = 2))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
#> [1] TRUE
```
