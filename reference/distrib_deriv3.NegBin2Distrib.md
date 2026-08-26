# Negative Binomial Third-Order Derivatives, NB2

Computes the four distinct third derivatives of the negative binomial
log-mass with respect to \\\mu\\ and \\\theta\\, in closed form. Writing
\\s = \theta + \mu\\, the three components involving the mean are
rational in \\(\mu, \theta)\\ and linear in \\y\\, and the pure
dispersion component carries \\\psi_2(y+\theta) - \psi_2(\theta)\\, with
\\\psi_2\\ the second derivative of the digamma function.

With `expected = TRUE` the expectations are returned. The three
components involving the mean need only \\\mathbb{E}\[Y\] = \mu\\; the
pure dispersion one needs \\\mathbb{E}\[\psi_2(Y+\theta)\]\\, which has
no closed form and is summed over the support with a far-tail
correction, exactly as the expected Hessian sums
\\\mathbb{E}\[\psi_1(Y+\theta)\]\\. `approx` and `nsim` are ignored
either way.

The cancellation of
[`distrib_gradient.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.NegBin2Distrib.md)
applies here too and **is not removed at this order**: the pure
dispersion components lose their digits at large \\\theta\\, where the
score itself is reliable to about five figures and no further.

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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_theta`,
`mu_theta_theta` and `theta_theta_theta`, each of length
`max(length(y), length(mu), length(theta))`. The names enumerate the
distinct multi-indices of order three in two parameters.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-mass in parameters
\\i\\, \\j\\ and \\k\\. \\\psi\\ is the digamma function and \\\psi_m\\
its \\m\\th derivative.

## See also

[`distrib_hessian.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin2Distrib.md)
for the order below and
[`distrib_deriv4.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.NegBin2Distrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- negbin2_distrib()
y <- c(0, 2, 6)
th <- list(mu = 4, theta = 2)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_theta"       "mu_theta_theta"   
#> [4] "theta_theta_theta"
d3$mu_mu_mu
#> [1] -0.01851852  0.02546296  0.11342593

# The expected values are constants at a fixed parameter setting.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] 0.06944444
#> 
#> $mu_mu_theta
#> [1] -0.02777778
#> 
#> $mu_theta_theta
#> [1] 0
#> 
#> $theta_theta_theta
#> [1] 0.08717095
#> 

# A central difference of the Hessian reproduces the same component.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 4 + eps, theta = 2))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 4 - eps, theta = 2))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
#> [1] TRUE
```
