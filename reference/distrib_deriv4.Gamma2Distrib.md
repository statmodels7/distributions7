# Gamma Fourth-Order Derivatives in Mean and Variance

Computes the five distinct fourth derivatives of the gamma log-density
with respect to \\\mu\\ and \\\sigma^2\\, in closed form, by the same
chain rule from the shape \\\alpha = \mu^2/\sigma^2\\ and the rate
\\\lambda = \mu/\sigma^2\\ that the lower orders use. Every component
carries \\\psi_3\\, the third derivative of the digamma function.

With `expected = TRUE` the expectations are returned, obtained by
replacing \\Y\\ with \\\mu\\ and \\\log Y\\ with \\\psi(\alpha) -
\log\lambda\\. Both routes are closed form, so `approx` and `nsim` are
ignored.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- y:

  A numeric vector of strictly positive observations. With
  `expected = TRUE` only its length is used.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

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

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_sigma2`,
`mu_mu_sigma2_sigma2`, `mu_sigma2_sigma2_sigma2` and
`sigma2_sigma2_sigma2_sigma2`, each of length
`max(length(y), length(mu), length(sigma2))`.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density in
parameters \\i\\, \\j\\, \\k\\ and \\l\\. \\\psi\\ is the digamma
function and \\\psi_m\\ its \\m\\th derivative.

## See also

[`distrib_deriv3.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gamma2Distrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic, and
[`distrib_expected_hessian.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma2Distrib.md)
for the second-order expectation these extend.

## Examples

``` r
d <- gamma2_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, sigma2 = 2)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"                 "mu_mu_mu_sigma2"            
#> [3] "mu_mu_sigma2_sigma2"         "mu_sigma2_sigma2_sigma2"    
#> [5] "sigma2_sigma2_sigma2_sigma2"

# The expected values are constants at a fixed parameter setting.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] 0.01065874
#> 
#> $mu_mu_mu_sigma2
#> [1] 0.1743748
#> 
#> $mu_mu_sigma2_sigma2
#> [1] -0.3074128
#> 
#> $mu_sigma2_sigma2_sigma2
#> [1] 0.259178
#> 
#> $sigma2_sigma2_sigma2_sigma2
#> [1] -0.5901462
#> 

# A central difference of the third order reproduces it.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(mu = 3 + eps, sigma2 = 2))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 3 - eps, sigma2 = 2))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-3)
#> [1] TRUE
```
