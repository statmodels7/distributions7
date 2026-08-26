# Gamma Third-Order Derivatives in Mean and Dispersion

Computes the four distinct third derivatives of the gamma log-density
with respect to \\\mu\\ and \\\phi\\, in closed form. With \\s =
1/\phi\\ and \\z = y/\mu\\, \$\$\ell^{(\mu\mu\mu)} = \dfrac{s(6z -
2)}{\mu^3}, \qquad \ell^{(\mu\mu\phi)} = -\dfrac{s^2(1 - 2z)}{\mu^2},
\qquad \ell^{(\mu\phi\phi)} = \dfrac{2s^3(z - 1)}{\mu},\$\$ and the pure
dispersion component follows from the derivatives in \\s\\ by Faa di
Bruno on \\s(\phi) = 1/\phi\\, \$\$\ell^{(\phi\phi\phi)} = f_3 (s')^3 +
3 f_2 s' s'' + f_1 s''',\$\$ with \\s' = -s^2\\, \\s'' = 2s^3\\, \\s'''
= -6s^4\\ and \\f_k\\ the \\k\\th derivative of \\\ell\\ in \\s\\.

With `expected = TRUE` the expectations are returned, obtained by
replacing \\z\\ with 1 and \\\log z\\ with \\\psi(s) - \log s\\: \\f_1\\
and the two components odd in \\z - 1\\ vanish. Both routes are closed
form, so no quadrature is run and `approx` and `nsim` are ignored.

## Arguments

- distrib:

  A `Gamma1Distrib` object, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

- y:

  A numeric vector of strictly positive observations. With
  `expected = TRUE` only its length is used.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  Both must be strictly positive.

- expected:

  Logical of length 1. When `TRUE` the expectations under the model are
  returned in place of the observed values. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available for both the observed and the
  expected values.

- nsim:

  Ignored, for the same reason. Defaults to `10000`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  the measured internal threshold the kernel stays sequential whatever
  the count says. Defaults to `1L`.

## Value

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_phi`,
`mu_phi_phi` and `phi_phi_phi`, each of length
`max(length(y), length(mu), length(phi))`. The names enumerate the
distinct multi-indices of order three in two parameters.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-density in
parameters \\i\\, \\j\\ and \\k\\. \\\psi\\ is the digamma function.

## See also

[`distrib_hessian.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Gamma1Distrib.md)
for the order below and
[`distrib_deriv4.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Gamma1Distrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic and for the numerical route a family without a closed
form takes.

## Examples

``` r
d <- gamma1_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, phi = 0.5)
d3 <- distrib_deriv3(d, y, th)

# The three components involving mu, written out.
s <- 1 / 0.5
z <- y / 3
all.equal(d3$mu_mu_mu, s * (6 * z - 2) / 3^3)
#> [1] TRUE
all.equal(d3$mu_mu_phi, -s^2 * (1 - 2 * z) / 3^2)
#> [1] TRUE
all.equal(d3$mu_phi_phi, 2 * s^3 * (z - 1) / 3)
#> [1] TRUE

# Expected values: the components odd in z - 1 vanish.
lapply(distrib_deriv3(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu
#> [1] 0.2962963
#> 
#> $mu_mu_phi
#> [1] 0.4444444
#> 
#> $mu_phi_phi
#> [1] 0
#> 
#> $phi_phi_phi
#> [1] 17.96406
#> 

# A central difference of the Hessian reproduces the same component.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 3, phi = 0.5 + eps))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 3, phi = 0.5 - eps))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_phi, tolerance = 1e-5)
#> [1] TRUE
```
