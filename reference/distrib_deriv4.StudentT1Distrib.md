# Student t Fourth-Order Derivatives

Computes the fifteen distinct fourth derivatives of the location-scale
Student t log-density in \\\mu\\, \\\sigma\\ and \\\nu\\. The observed
values are closed form and run in a compiled kernel decomposed over the
elements of the output, so they do not depend on the thread count.

**The expected values have no closed form.** With `expected = TRUE` the
method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
which integrates the observed derivatives against the density by the
strategy `approx` names. That is the one place on this page where
`approx` and `nsim` are read.

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
  returned in place of the value at the data, computed numerically.
  Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  One of `"integrate"` (the default here), `"bartlett"`, `"mc"` or
  `"opg"`, the strategy
  [`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
  uses. Read only when `expected = TRUE`.

- nsim:

  A single positive integer, the sample size when `approx = "mc"`. Read
  only when `expected = TRUE`. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Read
  only on the observed branch. Defaults to `1L`.

## Value

A named list of fifteen numeric vectors named for the multi-index they
carry, from `mu_mu_mu_mu` to `nu_nu_nu_nu`, each of length
`max(length(y), length(mu), length(sigma), length(nu))`.

## Large degrees of freedom

Unlike the third order, the fourth is **not** rewritten in the ratio
variables and ceases to be finite at a large \\\nu\\: measured at
\\\sigma = 1.2\\, eight of the fifteen components are finite at \\\nu =
10^{150}\\, five at \\10^{300}\\ and two at `.Machine$double.xmax`. A
`NaN` there is a loud failure and is preferable to a plausible wrong
number, and the regime is one in which the family is a Gaussian in all
but name. An outer criterion that reads this order at a \\\nu\\ run to
its clamp is reported as having no finite gradient rather than being
given one.

## See also

[`distrib_deriv3.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.StudentT1Distrib.md)
for the order below,
[`distrib_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentT1Distrib.md)
for the second order,
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
for the numerical expectation, and
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

# At a degrees of freedom the log link can produce, part of the order is
# not representable and says so.
big <- distrib_deriv4(d, y, list(mu = 0.4, sigma = 1.2, nu = 1e300))
sum(vapply(big, function(v) is.finite(v[1]), logical(1)))
#> [1] 5
```
