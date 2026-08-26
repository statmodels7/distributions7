# Skew Normal Third Derivatives

Computes the ten third derivatives of the log-density in closed form, in
a compiled kernel. The family's whole derivative surface stays
elementary for one reason: with \\t = \alpha z\\ and \\R(t) =
\phi(t)/\Phi(t)\\ the inverse Mills ratio, \\R' = -R(t+R)\\ closes the
recursion, so every derivative of \\\log\Phi(t)\\ is a polynomial in
\\t\\ and \\R\\ and no new special function appears at any order.

With `expected = TRUE` the value is an expectation instead, and there it
is **numerical**: the integrals are the ones that block the expected
information, so the call routes to
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
and `approx` and `nsim` are read.

## Arguments

- distrib:

  A `SkewNormal1Distrib` object, from
  [`skewnormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal1_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` its values
  are ignored and only its length matters, the result being one
  expectation repeated.

- theta:

  A named list with components `mu`, `sigma` and `alpha`, each a numeric
  vector of length 1 or of the length of `y`.

- expected:

  Logical of length 1. When `FALSE`, the default, the observed
  derivatives at `y` are returned from the compiled kernel. When `TRUE`
  the expectation is approximated by the strategy in `approx`.

- scale:

  Either `"parameter"`, the default, or `"link"`. The transformation is
  applied in the generic's body, so this method always returns the
  parameter scale.

- approx:

  One of `"integrate"`, `"bartlett"`, `"mc"` or `"opg"`, the strategy
  for the expectation. Read only when `expected = TRUE`; the first is
  the default and quadratures the observed kernel against the density.

- nsim:

  A single positive integer, the Monte Carlo sample size used when
  `approx = "mc"`. Ignored by every other strategy. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the compiled kernel may
  use. Defaults to `1L`. The result does not depend on the count.

## Value

A named list of ten numeric vectors, one per distinct third-order
component, named as
[`deriv_names()`](https://statmodels7.github.io/distributions7/reference/deriv_names.md)
names them: `mu_mu_mu`, `mu_mu_sigma`, `mu_mu_alpha`, `mu_sigma_sigma`,
`mu_sigma_alpha`, `mu_alpha_alpha`, `sigma_sigma_sigma`,
`sigma_sigma_alpha`, `sigma_alpha_alpha` and `alpha_alpha_alpha`.

## Notation

\\z = (y-\mu)/\sigma\\, \\t = \alpha z\\, \\R\\ the inverse Mills ratio
and \\\ell\\ the log-density of one observation.

## See also

[`distrib_hessian.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.SkewNormal1Distrib.md)
for the order below,
[`distrib_deriv4.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.SkewNormal1Distrib.md)
for the order above,
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- skewnormal1_distrib()
y <- c(-1.5, -0.3, 0.4, 2.1)
th <- list(mu = 0, sigma = 1, alpha = 3)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#>  [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_mu_alpha"      
#>  [4] "mu_sigma_sigma"    "mu_sigma_alpha"    "mu_alpha_alpha"   
#>  [7] "sigma_sigma_sigma" "sigma_sigma_alpha" "sigma_alpha_alpha"
#> [10] "alpha_alpha_alpha"

# Against a central difference of the analytic Hessian.
eps <- 1e-4
rbind(analytic = d3$mu_mu_alpha,
      numeric = (distrib_hessian(d, y, list(mu = 0, sigma = 1, alpha = 3 + eps))$mu_mu -
                 distrib_hessian(d, y, list(mu = 0, sigma = 1, alpha = 3 - eps))$mu_mu) /
                (2 * eps))
#>               [,1]      [,2]       [,3]         [,4]
#> analytic -5.953354 -5.070584 -0.8211369 6.658028e-07
#> numeric  -5.953354 -5.070584 -0.8211369 6.658030e-07

# At shape zero the third derivative in the location is the Gaussian's,
# which is exactly zero.
distrib_deriv3(d, y, list(mu = 0, sigma = 1, alpha = 0))$mu_mu_mu
#> [1] 0 0 0 0

# The expectation is numerical here, and the strategy is read.
set.seed(4)
c(integrate = distrib_deriv3(d, 0, th, expected = TRUE)$alpha_alpha_alpha,
  mc = distrib_deriv3(d, 0, th, expected = TRUE, approx = "mc",
                      nsim = 2000)$alpha_alpha_alpha)
#>  integrate         mc 
#> 0.01970008 0.02005633 
```
