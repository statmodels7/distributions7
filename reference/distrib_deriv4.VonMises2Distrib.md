# von Mises Fourth-Order Derivatives in the Resultant Length

Computes the five distinct fourth derivatives of the log-density in
\\\mu\\ and \\\rho\\, in closed form, by the construction
[`distrib_deriv3.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.VonMises2Distrib.md)
describes carried one order further: a single term \\D_a
\kappa^{(b)}(\rho)\\ for every component carrying a \\\mu\\, and the
fourth-order one-variable Faa di Bruno on \\\log I_0\\ for the
pure-\\\rho\\ one.

With `expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
which is the one place on this page where `approx` and `nsim` are read.

## Arguments

- distrib:

  A `VonMises2Distrib` object, from
  [`vonmises2_distrib()`](https://statmodels7.github.io/distributions7/reference/vonmises2_distrib.md).

- y:

  A numeric vector of angles in \\\[-\pi, \pi)\\. With `expected = TRUE`
  only its length is read.

- theta:

  A named list with components `mu` and `rho`, each a numeric vector of
  length 1 or of the length of `y`. `rho` must lie in \\(0, 1)\\.

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

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_rho`,
`mu_mu_rho_rho`, `mu_rho_rho_rho` and `rho_rho_rho_rho`, each of length
`length(y)`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\rho \in (0,1)\\ the mean resultant length, \\\kappa\\ the
concentration and \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\.

## See also

[`distrib_deriv3.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.VonMises2Distrib.md)
for the order below and the construction,
[`vm2_parts()`](https://statmodels7.github.io/distributions7/reference/vm2_parts.md)
for the map's derivatives, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d2 <- vonmises2_distrib()
y <- c(-1, 0, 0.5, 2)
th <- list(mu = 0.5, rho = 0.7)
d4 <- distrib_deriv4(d2, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"     "mu_mu_mu_rho"    "mu_mu_rho_rho"   "mu_rho_rho_rho" 
#> [5] "rho_rho_rho_rho"

# A central difference of the third order reproduces the pure-direction
# component.
eps <- 1e-5
up <- distrib_deriv3(d2, y, list(mu = 0.5 + eps, rho = 0.7))$mu_mu_mu
dn <- distrib_deriv3(d2, y, list(mu = 0.5 - eps, rho = 0.7))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-5)
#> [1] TRUE
```
