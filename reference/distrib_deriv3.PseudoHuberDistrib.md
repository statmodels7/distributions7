# Pseudo-Huber Third-Order Derivatives

Computes the ten distinct third derivatives of the pseudo-Huber
log-density in \\\mu\\, \\\sigma\\ and \\\nu\\, **in closed form**, in a
compiled kernel. Every component but the pure-\\\nu\\ one is a rational
function of \\r = y - \mu\\, \\\sigma\\ and \\D = \sqrt{\nu +
(r/\sigma)^2}\\; Bessel functions enter through \\\nu\\ alone, and the
exponentially scaled forms are used so that a large \\\nu\\ does not
overflow.

**The expected third derivatives have no closed form.** With
`expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
which integrates the observed derivatives against the density by the
strategy `approx` names. That is the one place on this page where
`approx` and `nsim` are read.

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

A named list of ten numeric vectors, `mu_mu_mu` through `nu_nu_nu`, each
of length `max(length(y), length(mu), length(sigma), length(nu))`.

## See also

[`distrib_hessian.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.PseudoHuberDistrib.md)
for the order below,
[`distrib_deriv4.PseudoHuberDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.PseudoHuberDistrib.md)
for the order above,
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
for the numerical expectation, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- pseudohuber_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 2)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#>  [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_mu_nu"         
#>  [4] "mu_sigma_sigma"    "mu_sigma_nu"       "mu_nu_nu"         
#>  [7] "sigma_sigma_sigma" "sigma_sigma_nu"    "sigma_nu_nu"      
#> [10] "nu_nu_nu"         

# A central difference of the Hessian reproduces the pure-location
# component, which is what says the closed form is the right one.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 0.4 + eps, sigma = 1.2, nu = 2))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, sigma = 1.2, nu = 2))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#> [1] TRUE

# The expected branch is a quadrature and takes a strategy; the components
# odd in the residual come back at zero.
vapply(distrib_deriv3(d, 0.4, th, expected = TRUE),
       function(v) v[1], numeric(1))
#>          mu_mu_mu       mu_mu_sigma          mu_mu_nu    mu_sigma_sigma 
#>     -6.938894e-18      3.055008e-01      2.612955e-02     -8.326673e-17 
#>       mu_sigma_nu          mu_nu_nu sigma_sigma_sigma    sigma_sigma_nu 
#>     -4.553649e-18     -1.192622e-18      3.135736e+00      7.878229e-02 
#>       sigma_nu_nu          nu_nu_nu 
#>      2.363469e-02      5.907273e-03 
```
