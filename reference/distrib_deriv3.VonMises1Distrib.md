# von Mises Third-Order Derivatives

Computes the four distinct third derivatives of the von Mises
log-density in \\\mu\\ and \\\kappa\\, in closed form. The log-density
is \\\kappa\cos(y-\mu) - \log\\2\pi I_0(\kappa)\\\\, **linear in**
\\\kappa\\ apart from the normalizing constant, so a component naming
one \\\mu\\ and two or more \\\kappa\\ is **exactly zero**. What remains
cycles: the pure-\\\mu\\ components run through \\\kappa\\\sin, -\cos,
-\sin, \cos\\(y-\mu)\\ with the order, the \\\mu\mu\kappa\\ component is
\\-\cos(y-\mu)\\, and the pure-\\\kappa\\ one is \\-A''(\kappa)\\, which
[`numericals7::bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.html)
supplies from the Riccati recursion \\A' = 1 - A/\kappa - A^2\\
differentiated.

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
  precedes `expected` here, unlike on most families, so both are best
  given by name.

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

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_kappa`,
`mu_kappa_kappa` and `kappa_kappa_kappa`, each of length `length(y)`.
`mu_kappa_kappa` is exactly zero.

## Notation

\\\ell\\ is the log-density of one observation, \\\mu\\ the mean
direction, \\\kappa \> 0\\ the concentration, \\I_m\\ the modified
Bessel function of the first kind of order \\m\\, and \\A(\kappa) =
I_1(\kappa)/I_0(\kappa)\\.

## See also

[`distrib_hessian.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.VonMises1Distrib.md)
for the order below,
[`distrib_deriv4.VonMises1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.VonMises1Distrib.md)
for the order above,
[`numericals7::bessel_i_ratio_derivs()`](https://statmodels7.github.io/numericals7/reference/bessel_i_ratio_derivs.html)
for the derivatives of \\A\\, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- vonmises1_distrib()
y <- c(-1, 0, 0.5, 2)
th <- list(mu = 0.5, kappa = 2)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "mu_mu_mu"          "mu_mu_kappa"       "mu_kappa_kappa"   
#> [4] "kappa_kappa_kappa"

# The component with one mu and two kappa is exactly zero, the log-density
# being linear in kappa apart from the constant.
d3$mu_kappa_kappa
#> [1] 0 0 0 0

# A central difference of the Hessian reproduces the pure-direction one.
eps <- 1e-5
up <- distrib_hessian(d, y, list(mu = 0.5 + eps, kappa = 2))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 0.5 - eps, kappa = 2))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-6)
#> [1] TRUE
```
