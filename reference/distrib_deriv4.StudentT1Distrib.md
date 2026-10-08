# Student t Fourth-Order Derivatives

Computes the fifteen distinct fourth derivatives of the location-scale
Student t log-density in \\\mu\\, \\\sigma\\ and \\\nu\\. The observed
values are closed form and run in a compiled kernel decomposed over the
elements of the output, so they do not depend on the thread count.

The expected values are closed forms, as for the third order (see
[`distrib_deriv3.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.StudentT1Distrib.md)).
`approx` and `nsim` are accepted for the generic's signature and
ignored.

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

A named list of fifteen numeric vectors named for the multi-index they
carry, from `mu_mu_mu_mu` to `nu_nu_nu_nu`, each of length
`max(length(y), length(mu), length(sigma), length(nu))`.

## Large degrees of freedom

Written like the third order, in \\z\\, \\q = z^2/\nu\\ and \\t =
1/(1+q)\\, with the quantities of \\\nu\\ alone on their asymptotic
series above \\\nu = 20\\: all fifteen components keep their digits and
stay finite to `.Machine$double.xmax`.

## See also

[`distrib_deriv3.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.StudentT1Distrib.md)
for the order below,
[`distrib_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentT1Distrib.md)
for the second order, and
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic.

## Examples

``` r
d <- student_t1_distrib()
y <- c(-2.5, 0.3, 1.8)
th <- list(mu = 0.4, sigma = 1.2, nu = 5)
d4 <- distrib_deriv4(d, y, th)
length(d4)
#> [1] 15
names(d4)[1:4]
#> [1] "mu_mu_mu_mu"       "mu_mu_mu_sigma"    "mu_mu_mu_nu"      
#> [4] "mu_mu_sigma_sigma"

# A central difference of the third order reproduces the pure-location
# component.
eps <- 1e-4
up <- distrib_deriv3(d, y, list(mu = 0.4 + eps, sigma = 1.2, nu = 5))$mu_mu_mu
dn <- distrib_deriv3(d, y, list(mu = 0.4 - eps, sigma = 1.2, nu = 5))$mu_mu_mu
all.equal((up - dn) / (2 * eps), d4$mu_mu_mu_mu, tolerance = 1e-5)
#> [1] TRUE

# At a degrees of freedom the log link can produce, every component is
# finite.
big <- distrib_deriv4(d, y, list(mu = 0.4, sigma = 1.2, nu = 1e300))
all(vapply(big, function(v) all(is.finite(v)), logical(1)))
#> [1] TRUE
```
