# Laplace Fifth-Order Derivatives

Computes the six distinct fifth derivatives of the Laplace log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form from a compiled
kernel. Away from the kink at \\y = \mu\\, with \\r = y - \mu\\, the
four components with two or more derivatives in \\\mu\\ are zero and
\$\$\ell^{(\mu\sigma\sigma\sigma\sigma)} =
\dfrac{24\\\mathrm{sign}(r)}{\sigma^5}, \qquad
\ell^{(\sigma\sigma\sigma\sigma\sigma)} = \dfrac{120 \|r\|}{\sigma^6} -
\dfrac{24}{\sigma^5}.\$\$ The kernel forms the powers of \\1/\sigma\\,
so that a small scale does not underflow an intermediate power where the
value is finite. With `scale = "link"` the components are composed with
the link derivatives by the routine that every family with a fifth-order
kernel uses.

## Arguments

- distrib:

  A `LaplaceDistrib` object, from
  [`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md).

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
d <- laplace_distrib()
y <- c(-1.2, 0.3, 2.5)
d5 <- distrib_deriv5(d, y, list(mu = 0.4, sigma = 1.5))
all.equal(d5$sigma_sigma_sigma_sigma_sigma,
          120 * abs(y - 0.4) / 1.5^6 - 24 / 1.5^5)
#> [1] TRUE
```
