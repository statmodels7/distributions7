# Pseudo-Huber Fourth-Order Derivatives

Computes the fifteen distinct fourth derivatives of the pseudo-Huber
log-density in \\\mu\\, \\\sigma\\ and \\\nu\\, in closed form, in a
compiled kernel. The Bessel functions enter through \\\nu\\ alone, and
their derivatives come from the three routes described in
[`distrib_deriv3.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.PseudoHuberDistrib.md).

The expected fourth derivatives have no closed form. With
`expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md).

## Arguments

- distrib:

  A `PseudoHuberDistrib` object, from
  [`pseudohuber_distrib()`](https://statmodels7.github.io/distributions7/reference/pseudohuber_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length is read.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectation under the model is
  returned in place of the value at the data, computed numerically.
  Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic reads it, and this method does not.

- method:

  `NULL` (the default) or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  It is read only when `expected` is `TRUE`: `NULL` selects
  [`expect_integrate()`](https://statmodels7.github.io/distributions7/reference/expect_integrate.md),
  and an object selects its own route. The generic applies an object
  before dispatching on the family, so that the object is used also by a
  family that has a closed form. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A named list of fifteen numeric vectors named for the multi-index that
they carry, from `mu_mu_mu_mu` to `nu_nu_nu_nu`, each of length
`max(length(y), length(mu), length(sigma), length(nu))`.

## See also

[`distrib_deriv3.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.PseudoHuberDistrib.md)
for the order below and the Bessel handling,
[`distrib_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md)
for the second order,
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
for the numerical expectation, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 2)
d4 <- distrib_deriv4(d, y, th)
length(d4)
#> [1] 15
names(d4)[1:4]
#> [1] "mu_mu_mu_mu"       "mu_mu_mu_sigma"    "mu_mu_mu_nu"      
#> [4] "mu_mu_sigma_sigma"

# A central difference of the third order reproduces the pure-location
# component.
eps <- 1e-4
up <- distrib_deriv3(d, y, list(mu = 0.4 + eps, sigma = 1.2, nu = 2))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 0.4 - eps, sigma = 1.2, nu = 2))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-5)
#> [1] TRUE
```
