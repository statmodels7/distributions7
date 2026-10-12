# Gaussian Fifth-Order Derivatives, Mean and Precision

Computes the six distinct fifth derivatives of the Gaussian log-density
with respect to \\\mu\\ and \\\tau\\, in closed form from a compiled
kernel. The log-density is quadratic in \\\mu\\ and linear in \\\tau\\
apart from \\\tfrac12\log\tau\\, so every component is zero except
\$\$\ell^{(\tau\tau\tau\tau\tau)} = \dfrac{12}{\tau^5}.\$\$ With
`scale = "link"` the components are composed with the link derivatives
of both parameters by the same routine that the families with a
fifth-order kernel use.

## Arguments

- distrib:

  A `Gaussian3Distrib` object, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `tau`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.

- scale:

  One of `"parameter"` (the default) or `"link"`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use.
  Defaults to `1L`.

## Value

A named list of six numeric vectors, `mu_mu_mu_mu_mu`,
`mu_mu_mu_mu_tau`, `mu_mu_mu_tau_tau`, `mu_mu_tau_tau_tau`,
`mu_tau_tau_tau_tau` and `tau_tau_tau_tau_tau`, each of length
`max(length(y), length(mu), length(tau))`.

## See also

[`distrib_deriv4.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian3Distrib.md)
for the order below and
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
for the generic.

## Examples

``` r
d <- gaussian3_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, tau = 1.5)
d5 <- distrib_deriv5(d, y, th)
names(d5)
#> [1] "mu_mu_mu_mu_mu"      "mu_mu_mu_mu_tau"     "mu_mu_mu_tau_tau"   
#> [4] "mu_mu_tau_tau_tau"   "mu_tau_tau_tau_tau"  "tau_tau_tau_tau_tau"
all.equal(d5$tau_tau_tau_tau_tau, rep(12 / 1.5^5, 3))
#> [1] TRUE
```
