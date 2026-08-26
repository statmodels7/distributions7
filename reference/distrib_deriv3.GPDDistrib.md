# Generalized Pareto Third-Order Derivatives

Computes the four distinct third derivatives of the generalized Pareto
log-density in \\\sigma\\ and \\\xi\\, **in closed form**, through
[`gpd_components()`](https://statmodels7.github.io/distributions7/reference/gpd_components.md).
The log-density splits as \\-\log\sigma - \log t - \log(t)/\xi\\ with
\\t = 1 + \xi y/\sigma\\, and the last piece is taken from its own
series wherever the Leibniz form's terms of size \\\xi^{-(b+1)}\\ cancel
against each other.

With `expected = TRUE` the method calls
[`expected_derivative()`](https://statmodels7.github.io/distributions7/reference/expected_derivative.md)
instead: the expected third derivatives have no closed form. That is the
one place on this page where `approx` and `nsim` are read.

## Arguments

- distrib:

  A `GPDDistrib` object, from
  [`gpd_distrib()`](https://statmodels7.github.io/distributions7/reference/gpd_distrib.md).

- y:

  A numeric vector of observations on the support. With
  `expected = TRUE` only its length is read.

- theta:

  A named list with components `sigma` and `xi`, each a numeric vector
  of length 1 or of the length of `y`. `sigma` must be strictly
  positive; `xi` may be of either sign, including zero.

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

- threads:

  A single positive integer, how many threads the polynomial kernel of
  the series branch may use. Defaults to `1L`.

## Value

A named list of four numeric vectors, `sigma_sigma_sigma`,
`sigma_sigma_xi`, `sigma_xi_xi` and `xi_xi_xi`, each of length
`max(length(y), lengths(theta))`.

## Notation

\\\ell\\ is the log-density of one observation, \\\sigma \> 0\\ the
scale, \\\xi\\ the shape, \\z = y/\sigma\\ and \\t = 1 + \xi z\\.

## See also

[`distrib_hessian.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.GPDDistrib.md)
for the order below,
[`distrib_deriv4.GPDDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.GPDDistrib.md)
for the order above,
[`gpd_components()`](https://statmodels7.github.io/distributions7/reference/gpd_components.md)
for the two routes and the measured threshold, and
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- gpd_distrib()
y <- c(0.5, 2, 8)
th <- list(sigma = 1.5, xi = 0.3)
d3 <- distrib_deriv3(d, y, th)
names(d3)
#> [1] "sigma_sigma_sigma" "sigma_sigma_xi"    "sigma_xi_xi"      
#> [4] "xi_xi_xi"         

# A central difference of the Hessian reproduces the pure-scale component.
eps <- 1e-5
up <- distrib_hessian(d, y, list(sigma = 1.5 + eps, xi = 0.3))$sigma_sigma
dn <- distrib_hessian(d, y, list(sigma = 1.5 - eps, xi = 0.3))$sigma_sigma
all.equal((up - dn) / (2 * eps), d3$sigma_sigma_sigma, tolerance = 1e-6)
#> [1] TRUE

# At a shape near zero the family is the exponential, and the scale
# components converge onto that family's at rate O(xi).
de <- exponential_distrib()
vapply(c(1e-6, 1e-9), function(x)
  max(abs(distrib_deriv3(d, y, list(sigma = 1.5, xi = x))$sigma_sigma_sigma -
          distrib_deriv3(de, y, list(mu = 1.5))$mu_mu_mu)), numeric(1))
#> [1] 9.165352e-05 9.165432e-08
```
