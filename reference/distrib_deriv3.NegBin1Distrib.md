# NB1 Third-Order Derivatives

Computes the four distinct third derivatives of the NB1 log-likelihood
in the mean \\\mu\\ and the dispersion \\\theta\\, **in closed form**,
through
[`negbin1_components()`](https://statmodels7.github.io/distributions7/reference/negbin1_components.md).
In the size \\r = \mu/\theta\\ the log-likelihood is \\G(r) +
rB(\theta) + C(\theta)\\, so the only composite piece is
\\G(\mu/\theta)\\, and its mixed derivatives follow a recursion in the
powers of \\r\\ and the order of \\G\\ that is run rather than solved.

With `expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
instead: the expected third derivatives have no closed form. That is the
one place on this page where `approx` and `nsim` are read.

## Arguments

- distrib:

  A `NegBin1Distrib` object, from
  [`negbin1_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin1_distrib.md).

- y:

  A numeric vector of counts. With `expected = TRUE` only its length is
  read.

- theta:

  A named list with components `mu` and `theta`, each a numeric vector
  of length 1 or of the length of `y`, both strictly positive. A
  component of length 1 is recycled.

- expected:

  Logical of length 1. When `TRUE` the expectation under the model is
  returned in place of the value at the data, computed numerically.
  Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  One of `"integrate"` (the default here), `"bartlett"`, `"mc"` or
  `"opg"`. Read only when `expected = TRUE`.

- nsim:

  A single positive integer, the sample size when `approx = "mc"`. Read
  only when `expected = TRUE`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_theta`,
`mu_theta_theta` and `theta_theta_theta`, each of length
`max(length(y), lengths(theta))`.

## Notation

\\\ell\\ is the log-likelihood of one observation, \\\mu \> 0\\ the
mean, \\\theta \> 0\\ the dispersion, \\r = \mu/\theta\\ the negative
binomial size and \\G(r) = \log\Gamma(y+r) - \log\Gamma(r)\\.

## The Poisson boundary

As \\\theta \to 0\\ the family tends to the Poisson and the recursion's
terms in the powers of \\r\\ grow while their sum stays of order one.
The polygamma differences go through
[`psi_shift_diff()`](https://statmodels7.github.io/distributions7/reference/psi_shift_diff.md)
and are exact, but the cancellation among those powers is not repaired:
at \\\theta = 5\times10^{-4}\\ the terms reach \\8\times10^{6}\\, and
this order is not reliable below about \\\theta = 0.05\\. The score
itself, which does not carry the recursion, reaches the Poisson limit to
five figures.

## See also

[`distrib_hessian.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.NegBin1Distrib.md)
for the order below,
[`distrib_deriv4.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.NegBin1Distrib.md)
for the order above,
[`negbin1_components()`](https://statmodels7.github.io/distributions7/reference/negbin1_components.md)
for the recursion,
[`distrib_deriv3.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin2Distrib.md)
for the other negative binomial, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- negbin1_distrib()
y <- c(0, 3, 7)
th <- list(mu = 4, theta = 1.2)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_theta"       "mu_theta_theta"   
#> [4] "theta_theta_theta"

# A central difference of the Hessian reproduces the pure-mean component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 4 + eps, theta = 1.2))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 4 - eps, theta = 1.2))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#> [1] TRUE

# And a mixed component, which is where the recursion does its work.
up <- distrib_hessian(d, y, list(mu = 4, theta = 1.2 + eps))$mu_theta
dn <- distrib_hessian(d, y, list(mu = 4, theta = 1.2 - eps))$mu_theta
all.equal((up - dn) / (2 * eps), d3$mu_theta_theta, tolerance = 1e-6)
#> [1] TRUE
```
