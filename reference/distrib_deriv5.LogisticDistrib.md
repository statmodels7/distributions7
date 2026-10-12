# Logistic Fifth-Order Derivatives

Computes the six distinct fifth derivatives of the logistic log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form from a compiled
kernel. Writing the log-density as \\\ell = -\log\sigma + g(z)\\ with
\\z = (y - \mu)/\sigma\\, each component is \\\sigma^{-5}\\ times a sum
of terms \\c\\ z^p g^{(i)}(z)\\ with integer coefficients, plus
\\-24\sigma^{-5}\\ in the component that is pure in \\\sigma\\. The
kernel evaluates the sum from the derivatives of \\g(z) = -z - 2\log(1 +
e^{-z})\\, written in the sigmoid \\t = 1/(1 + e^{-z})\\ and \\u = 1 -
t\\, so that every factor \\g^{(i)}\\ lies in \\\[-1, 1\]\\. With
`scale = "link"` the components are composed with the link derivatives
by the routine that every family with a fifth-order kernel uses.

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`.

- scale:

  One of `"parameter"` (the default) or `"link"`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of six numeric vectors, `mu_mu_mu_mu_mu`,
`mu_mu_mu_mu_sigma`, `mu_mu_mu_sigma_sigma`, `mu_mu_sigma_sigma_sigma`,
`mu_sigma_sigma_sigma_sigma` and `sigma_sigma_sigma_sigma_sigma`, each
of length `max(length(y), length(mu), length(sigma))`.

## See also

[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the order below and
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
d5 <- distrib_deriv5(d, c(-1.2, 0.3, 2.5), list(mu = 0.4, sigma = 1.5))
names(d5)
#> [1] "mu_mu_mu_mu_mu"                "mu_mu_mu_mu_sigma"            
#> [3] "mu_mu_mu_sigma_sigma"          "mu_mu_sigma_sigma_sigma"      
#> [5] "mu_sigma_sigma_sigma_sigma"    "sigma_sigma_sigma_sigma_sigma"
```
