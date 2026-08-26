# Negative Binomial Fourth-Order Derivatives, NB2

Computes the five distinct fourth derivatives of the negative binomial
log-mass with respect to \\\mu\\ and \\\theta\\, in closed form, by the
same route the third order takes: the components involving the mean are
rational in \\(\mu, \theta)\\ and linear in \\y\\, and the pure
dispersion component carries \\\psi_3(y+\theta) - \psi_3(\theta)\\, with
\\\psi_3\\ the third derivative of the digamma function.

With `expected = TRUE` the expectations are returned. The pure
dispersion one needs \\\mathbb{E}\[\psi_3(Y+\theta)\]\\, summed over the
support with a far-tail correction; the rest need only \\\mathbb{E}\[Y\]
= \mu\\. `approx` and `nsim` are ignored either way, and the caveat
about large \\\theta\\ on
[`distrib_deriv3.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin2Distrib.md)
applies here as well.

## Arguments

- distrib:

  A `NegBin2Distrib` object, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- y:

  A numeric vector of counts. With `expected = TRUE` only its length is
  used.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
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

  Ignored, the expected values being closed form or exact sums.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_theta`,
`mu_mu_theta_theta`, `mu_theta_theta_theta` and
`theta_theta_theta_theta`, each of length
`max(length(y), length(mu), length(theta))`.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-mass in parameters
\\i\\, \\j\\, \\k\\ and \\l\\. \\\psi\\ is the digamma function and
\\\psi_m\\ its \\m\\th derivative.

## See also

[`distrib_deriv3.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin2Distrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic, and
[`distrib_expected_hessian.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.NegBin2Distrib.md)
for the second-order expectation these extend.

## Examples

``` r
d <- negbin2_distrib()
y <- c(0, 2, 6)
th <- list(mu = 4, theta = 2)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_theta"         
#> [3] "mu_mu_theta_theta"       "mu_theta_theta_theta"   
#> [5] "theta_theta_theta_theta"

# The expected values are constants at a fixed parameter setting.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] -0.06597222
#> 
#> $mu_mu_mu_theta
#> [1] 0.01851852
#> 
#> $mu_mu_theta_theta
#> [1] 0.009259259
#> 
#> $mu_theta_theta_theta
#> [1] 0
#> 
#> $theta_theta_theta_theta
#> [1] -0.1686028
#> 

# A central difference of the third order reproduces it.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(mu = 4 + eps, theta = 2))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 4 - eps, theta = 2))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-3)
#> [1] TRUE
```
