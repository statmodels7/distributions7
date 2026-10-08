# Student t Third-Order Derivatives

Computes the ten distinct third derivatives of the location-scale
Student t log-density in \\\mu\\, \\\sigma\\ and \\\nu\\. The observed
values are closed form and run in a compiled kernel decomposed over the
elements of the output, so they do not depend on the thread count.

The expected values are closed forms too: under the model \\1/(1 +
z^2/\nu)\\ follows a beta distribution with parameters \\\nu/2\\ and
\\1/2\\, so the expectation of every component is a rational function of
\\\nu\\ plus the derivatives of \\\log\Gamma((\nu+1)/2) -
\log\Gamma(\nu/2)\\. `approx` and `nsim` are accepted for the generic's
signature and ignored.

## Arguments

- distrib:

  A `StudentT1Distrib` object, from
  [`student_t1_distrib()`](https://statmodels7.github.io/distributions7/reference/student_t1_distrib.md).

- y:

  A numeric vector of observations. With `expected = TRUE` only its
  length is read.

- theta:

  A named list with components `mu`, `sigma` and `nu`, each a numeric
  vector of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` and `nu` must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectation under the model is
  returned in place of the value at the data. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx, nsim:

  Accepted for the generic's signature; the expectations are closed
  forms.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of ten numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
`mu_mu_nu`, `mu_sigma_sigma`, `mu_sigma_nu`, `mu_nu_nu`,
`sigma_sigma_sigma`, `sigma_sigma_nu`, `sigma_nu_nu` and `nu_nu_nu`,
each of length `max(length(y), length(mu), length(sigma), length(nu))`.

## Large degrees of freedom

See
[`distrib_gradient.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.StudentT1Distrib.md):
every order is written in \\z = r/\sigma\\, \\q = z^2/\nu\\ and \\t =
1/(1+q)\\, and the quantities of \\\nu\\ alone switch to their
asymptotic series above \\\nu = 20\\, so all ten components keep their
digits and stay finite to `.Machine$double.xmax`. On the link scale the
chain rule forms \\(h')^k\\ against a component of order \\\nu^{-k}\\,
and one of the ten ceases to be finite at \\\nu = 10^{150}\\.

## See also

[`distrib_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentT1Distrib.md)
for the order below,
[`distrib_deriv4.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.StudentT1Distrib.md)
for the order above, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 5)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#>  [1] "mu_mu_mu"          "mu_mu_sigma"       "mu_mu_nu"         
#>  [4] "mu_sigma_sigma"    "mu_sigma_nu"       "mu_nu_nu"         
#>  [7] "sigma_sigma_sigma" "sigma_sigma_nu"    "sigma_nu_nu"      
#> [10] "nu_nu_nu"         

# A central difference of the Hessian reproduces the pure-location
# component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 0.4 + eps, sigma = 1.2, nu = 5))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, sigma = 1.2, nu = 5))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#> [1] TRUE

# Averaging the observed branch over draws reaches the expected one; the
# components odd in the residual go to zero.
set.seed(2)
z <- distrib_rng(d, 2e5, th)
rbind(expected = vapply(distrib_deriv3(d, y, th, expected = TRUE),
                        function(v) v[1], numeric(1)),
      averaged = vapply(distrib_deriv3(d, z, th), mean, numeric(1)))
#>              mu_mu_mu mu_mu_sigma    mu_mu_nu mu_sigma_sigma  mu_sigma_nu
#> expected 0.0000000000   0.6076389 -0.01388889    0.000000000 0.0000000000
#> averaged 0.0007771678   0.6061054 -0.01396171   -0.003885839 0.0001479138
#>              mu_nu_nu sigma_sigma_sigma sigma_sigma_nu  sigma_nu_nu    nu_nu_nu
#> expected 0.000000e+00          2.748843     -0.1041667 -0.005555556 0.002511281
#> averaged 8.536436e-06          2.756510     -0.1044698 -0.005578214 0.002510767
```
