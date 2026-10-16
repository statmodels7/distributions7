# Cauchy Fifth-Order Derivatives

Computes the six distinct fifth derivatives of the Cauchy log-density
with respect to \\\mu\\ and \\\sigma\\, in closed form from a compiled
kernel. With \\Y = y - \mu + \mathrm{i}\sigma\\ the log-density is
\\\ell = \log\sigma - 2\\\mathrm{Re}\log Y - \log\pi\\, and the
derivative of \\\log Y\\ of order \\n\\, with \\j\\ derivatives in
\\\mu\\ and \\k = n - j\\ in \\\sigma\\, is \\(-1)^j \mathrm{i}^k
(-1)^{n-1}(n-1)!\\Y^{-n}\\. Each component is therefore the single term
\$\$-2(-1)^{j+n-1}(n-1)!\\\mathrm{Re}(\mathrm{i}^k Y^{-n}),\$\$ plus
\\(-1)^{n-1}(n-1)!/\sigma^n\\ from \\\log\sigma\\ in the component that
is pure in \\\sigma\\. The kernel forms \\1/Y\\ by Smith's division and
never forms \\(y - \mu)^2\\ or \\\sigma^2\\, so the terms do not cancel
one another and the components stay finite at any residual. The observed
orders one to four are computed by the same formula. With
`scale = "link"` the components are composed with the link derivatives
by the routine that every family with a fifth-order kernel uses.

## Arguments

- distrib:

  A `CauchyDistrib` object, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

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
d <- cauchy_distrib()
d5 <- distrib_deriv5(d, c(-1.2, 0.3, 1e200), list(mu = 0.4, sigma = 1.5))
# The derivatives are finite far in the tail, where the component pure in
# sigma tends to 24/sigma^5.
c(d5$sigma_sigma_sigma_sigma_sigma[3], 24 / 1.5^5)
#> [1] 3.160494 3.160494
```
