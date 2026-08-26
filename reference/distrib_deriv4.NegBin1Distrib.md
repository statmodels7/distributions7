# NB1 Fourth-Order Derivatives

Computes the five distinct fourth derivatives of the NB1 log-likelihood
in \\\mu\\ and \\\theta\\, **in closed form**, by the construction
[`distrib_deriv3.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin1Distrib.md)
describes carried one order further: the same coefficient recursion over
the powers of \\r = \mu/\theta\\ and the order of \\G(r) =
\log\Gamma(y+r) - \log\Gamma(r)\\.

With `expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
instead: the expected fourth derivatives have no closed form.

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

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_theta`,
`mu_mu_theta_theta`, `mu_theta_theta_theta` and
`theta_theta_theta_theta`, each of length
`max(length(y), lengths(theta))`.

## Notation

\\\ell\\ is the log-likelihood of one observation, \\\mu \> 0\\ the
mean, \\\theta \> 0\\ the dispersion, \\r = \mu/\theta\\ the negative
binomial size and \\G(r) = \log\Gamma(y+r) - \log\Gamma(r)\\.

## The Poisson boundary

The caveat of
[`distrib_deriv3.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin1Distrib.md)
applies here and more strongly: the recursion divides by \\\theta^{4}\\,
so the cancellation among the powers of \\r\\ is worse at this order
than at the one below. This order is not reliable below about \\\theta =
0.05\\.

## See also

[`distrib_deriv3.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.NegBin1Distrib.md)
for the order below and the recursion,
[`negbin1_components()`](https://statmodels7.github.io/distributions7/reference/negbin1_components.md)
for the assembly, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- negbin1_distrib()
y <- c(0, 3, 7)
th <- list(mu = 4, theta = 1.2)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_theta"         
#> [3] "mu_mu_theta_theta"       "mu_theta_theta_theta"   
#> [5] "theta_theta_theta_theta"

# A central difference of the third order reproduces a mixed component.
eps <- 1e-4
up <- distrib_deriv3(d, y, list(mu = 4, theta = 1.2 + eps))$mu_mu_theta
dn <- distrib_deriv3(d, y, list(mu = 4, theta = 1.2 - eps))$mu_mu_theta
all.equal((up - dn) / (2 * eps), d4$mu_mu_theta_theta, tolerance = 1e-5)
#> [1] TRUE
```
