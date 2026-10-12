# Gaussian Fifth-Order Derivatives

Computes the six distinct fifth derivatives of the Gaussian log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form from a compiled
kernel. Writing \\r = y - \mu\\, the three components with three or more
derivatives in \\\mu\\ are zero and
\$\$\ell^{(\mu\mu\sigma\sigma\sigma)} = \dfrac{24}{\sigma^5}, \qquad
\ell^{(\mu\sigma\sigma\sigma\sigma)} = \dfrac{120 r}{\sigma^6}, \qquad
\ell^{(\sigma\sigma\sigma\sigma\sigma)} = -\dfrac{24}{\sigma^5} +
\dfrac{360 r^2}{\sigma^7}.\$\$ With `scale = "link"` the components are
composed with the link derivatives of both parameters by the same
routine that the families with a fifth-order kernel use.

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

- y:

  A numeric vector of observations.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled.

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

[`distrib_deriv4.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gaussian1Distrib.md)
for the order below and
[`distrib_deriv5()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv5.md)
for the generic.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1.2, 0.3, 2.5)
th <- list(mu = 0.4, sigma = 1.5)
d5 <- distrib_deriv5(d, y, th)
names(d5)
#> [1] "mu_mu_mu_mu_mu"                "mu_mu_mu_mu_sigma"            
#> [3] "mu_mu_mu_sigma_sigma"          "mu_mu_sigma_sigma_sigma"      
#> [5] "mu_sigma_sigma_sigma_sigma"    "sigma_sigma_sigma_sigma_sigma"
all.equal(d5$sigma_sigma_sigma_sigma_sigma,
          -24 / 1.5^5 + 360 * (y - 0.4)^2 / 1.5^7)
#> [1] TRUE
```
