# von Mises Fourth-Order Derivatives

Computes the five distinct fourth derivatives of the von Mises
log-density in \\\mu\\ and \\\kappa\\, in closed form, by the
construction
[`distrib_deriv3.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.VonMises1Distrib.md)
describes carried one order further. The log-density is linear in
\\\kappa\\ apart from the normalizing constant, so **two of the five are
exactly zero**: `mu_mu_kappa_kappa` and `mu_kappa_kappa_kappa`. The
pure-direction component is \\\kappa\cos(y-\mu)\\, the
\\\mu\mu\mu\kappa\\ one \\-\sin(y-\mu)\\, and the pure-concentration one
\\-A^{(3)}(\kappa)\\.

With `expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
which is the one place on this page where `approx` and `nsim` are read.

## Arguments

- distrib:

  A `VonMises1Distrib` object, from
  [`vonmises1_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises1_distrib.md).

- y:

  A numeric vector of angles in \\\[-\pi, \pi)\\. With `expected = TRUE`
  only its length is read.

- theta:

  A named list with components `mu` and `kappa`, each a numeric vector
  of length 1 or of the length of `y`. `kappa` must be strictly
  positive.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method. **Note the argument order**: `scale`
  precedes `expected` here, so both are best given by name.

- expected:

  Logical of length 1. When `TRUE` the expectation under the model is
  returned in place of the value at the data, computed numerically.
  Defaults to `FALSE`.

- approx:

  One of `"integrate"` (the default here), `"bartlett"`, `"mc"` or
  `"opg"`. Read only when `expected = TRUE`.

- nsim:

  A single positive integer, the sample size when `approx = "mc"`. Read
  only when `expected = TRUE`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_kappa`,
`mu_mu_kappa_kappa`, `mu_kappa_kappa_kappa` and
`kappa_kappa_kappa_kappa`, each of length `length(y)`. The middle two
are exactly zero.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\kappa \> 0\\ the concentration, and \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\.

## See also

[`distrib_deriv3.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.VonMises1Distrib.md)
for the order below and the construction,
[`numericals7::bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.html)
for the derivatives of \\A\\, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- vonmises1_distrib()
y <- c(-1, 0, 0.5, 2)
th <- list(mu = 0.5, kappa = 2)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"             "mu_mu_mu_kappa"         
#> [3] "mu_mu_kappa_kappa"       "mu_kappa_kappa_kappa"   
#> [5] "kappa_kappa_kappa_kappa"

# Two of the five are exactly zero.
c(d4$mu_mu_kappa_kappa[1], d4$mu_kappa_kappa_kappa[1])
#> [1] 0 0

# A central difference of the third order reproduces the pure-direction
# component.
eps <- 1e-4
up <- distrib_deriv3(d, y, list(mu = 0.5 + eps, kappa = 2))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 0.5 - eps, kappa = 2))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-5)
#> [1] TRUE
```
