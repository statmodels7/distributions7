# Beta Fourth-Order Derivatives in Mean and Precision

Computes the five distinct fourth derivatives of the beta log-density
with respect to \\\mu\\ and \\\phi\\, in closed form as combinations of
\\\psi_3\\, the third derivative of the digamma function, at \\\alpha =
\mu\phi\\, \\\beta = (1-\mu)\phi\\ and \\\phi\\.

As at third order, every component is free of the response, so the
observed and expected values coincide exactly and `expected`, `approx`
and `nsim` are all without effect.

## Arguments

- distrib:

  A `Beta1Distrib` object, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `mu` must lie strictly in \\(0, 1)\\ and `phi` must be strictly
  positive.

- expected:

  Logical of length 1, and without effect here, the observed and
  expected fourth derivatives being the same numbers. Defaults to
  `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available.

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
`length(y)` and constant within itself when the parameters are.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density in
parameters \\i\\, \\j\\, \\k\\ and \\l\\. \\\psi\\ is the digamma
function and \\\psi_m\\ its \\m\\th derivative.

## See also

[`distrib_deriv3.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Beta1Distrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic, and
[`distrib_expected_hessian.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Beta1Distrib.md)
for the second-order expectation these extend.

## Examples

``` r
d <- beta1_distrib()
y <- c(0.2, 0.5, 0.8)
th <- list(mu = 0.4, phi = 5)
d4 <- distrib_deriv4(d, y, th)

# Five constants: nothing at this order depends on the observation.
lapply(d4, unique)
#> $mu_mu_mu_mu
#> [1] -383.0493
#> 
#> $mu_mu_mu_phi
#> [1] 2.973485
#> 
#> $mu_mu_phi_phi
#> [1] -0.04367237
#> 
#> $mu_phi_phi_phi
#> [1] -0.002074338
#> 
#> $phi_phi_phi_phi
#> [1] -0.006631567
#> 

# So asking for the expectation changes nothing.
identical(d4, distrib_deriv4(d, y, th, expected = TRUE))
#> [1] TRUE

# A central difference of the third order reproduces it.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(mu = 0.4 + eps, phi = 5))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 0.4 - eps, phi = 5))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-3)
#> [1] TRUE
```
