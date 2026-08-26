# von Mises Third-Order Derivatives in the Resultant Length

Computes the four distinct third derivatives of the log-density in
\\\mu\\ and \\\rho\\, in closed form. Every component carrying at least
one \\\mu\\ collapses to a **single term** \\D_a \kappa^{(b)}(\rho)\\,
with \\D_a\\ the \\a\\-th \\\mu\\-derivative of \\\cos(y-\mu)\\ and
\\\kappa^{(b)}\\ the \\b\\-th derivative of \\A^{-1}\\: the
concentration parametrization's \\\mu\\-derivatives are linear in
\\\kappa\\, so the composition has nothing to expand. The pure-\\\rho\\
component carries the full one-variable Faa di Bruno on \\\log I_0\\,
written out.

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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_rho`,
`mu_rho_rho` and `rho_rho_rho`, each of length `length(y)`.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\rho \in (0,1)\\ the mean resultant length, \\\kappa\\ the
concentration and \\A(\kappa) = I_1(\kappa)/I_0(\kappa)\\.

## See also

[`distrib_hessian.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises2Distrib.md)
for the order below,
[`distrib_deriv4.VonMises2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.VonMises2Distrib.md)
for the order above,
[`vm2_parts()`](https://statmodels7.github.io/distributions7/reference/vm2_parts.md)
for the map's derivatives, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d2 <- vonmises2_distrib()
y <- c(-1, 0, 0.5, 2)
th <- list(mu = 0.5, rho = 0.7)
d3 <- distrib_deriv3(d2, y, th)
names(d3)
#> [1] "mu_mu_mu"    "mu_mu_rho"   "mu_rho_rho"  "rho_rho_rho"

# A central difference of the Hessian reproduces the pure-rho component,
# which is the one carrying the whole chain rule.
eps <- 1e-5
up <- distrib_hessian(d2, y, list(mu = 0.5, rho = 0.7 + eps))$rho_rho
dn <- distrib_hessian(d2, y, list(mu = 0.5, rho = 0.7 - eps))$rho_rho
all.equal((up - dn) / (2 * eps), d3$rho_rho_rho, tolerance = 1e-5)
#> [1] TRUE

# Unlike the concentration parametrization, no component is exactly zero:
# the map's higher derivatives bring the data into all four.
vapply(d3, function(v) v[1], numeric(1))
#>     mu_mu_mu    mu_mu_rho   mu_rho_rho  rho_rho_rho 
#>    2.0085836   -0.4356578  -31.5311208 -267.6200994 
```
