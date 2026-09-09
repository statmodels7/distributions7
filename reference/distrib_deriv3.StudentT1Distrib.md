# Student t Third-Order Derivatives

Computes the ten distinct third derivatives of the location-scale
Student t log-density in \\\mu\\, \\\sigma\\ and \\\nu\\. The observed
values are closed form and run in a compiled kernel decomposed over the
elements of the output, so they do not depend on the thread count.

**The expected values have no closed form.** With `expected = TRUE` the
method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md),
which integrates the observed derivatives against the density by the
strategy `approx` names. That is the one place on this page where
`approx` and `nsim` are read; on the observed branch both are ignored.

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

A named list of ten numeric vectors, `mu_mu_mu`, `mu_mu_sigma`,
`mu_mu_nu`, `mu_sigma_sigma`, `mu_sigma_nu`, `mu_nu_nu`,
`sigma_sigma_sigma`, `sigma_sigma_nu`, `sigma_nu_nu` and `nu_nu_nu`,
each of length `max(length(y), length(mu), length(sigma), length(nu))`.

## Large degrees of freedom

Every component is divided by \\D^3\\ with \\D = \nu\sigma^2 + r^2\\,
and \\D^3\\ overflows at \\5.6\times10^{102}\\ where the log link
reaches \\1.8\times10^{308}\\. The shipped kernel is written in \\z =
r/\sigma\\, \\u = z^2/\nu\\ and \\t = 1/(1+u)\\ instead, so all ten stay
finite to `.Machine$double.xmax`. **On the link scale they do not**: the
chain rule forms \\(h')^k\\ against a component of order \\\nu^{-k}\\,
and one of the ten ceases to be finite at \\\nu = 10^{150}\\. That
regime is where the family is a Gaussian in all but name.

## See also

[`distrib_hessian.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.StudentT1Distrib.md)
for the order below,
[`distrib_deriv4.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.StudentT1Distrib.md)
for the order above,
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
for the numerical expectation, and
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

# The expected branch is a quadrature, and averaging the observed one over
# draws reaches it; the components odd in the residual go to zero.
set.seed(2)
z <- distrib_rng(d, 2e5, th)
rbind(expected = vapply(distrib_deriv3(d, y, th, expected = TRUE),
                        function(v) v[1], numeric(1)),
      averaged = vapply(distrib_deriv3(d, z, th), mean, numeric(1)))
#>              mu_mu_mu mu_mu_sigma    mu_mu_nu mu_sigma_sigma   mu_sigma_nu
#> expected 2.775558e-17   0.6076389 -0.01388889  -1.942890e-16 -3.469447e-18
#> averaged 7.771678e-04   0.6061054 -0.01396171  -3.885839e-03  1.479138e-04
#>              mu_nu_nu sigma_sigma_sigma sigma_sigma_nu  sigma_nu_nu    nu_nu_nu
#> expected 2.168404e-19          2.748843     -0.1041667 -0.005555556 0.002511281
#> averaged 8.536436e-06          2.756510     -0.1044698 -0.005578214 0.002510767
```
