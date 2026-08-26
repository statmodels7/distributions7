# Gamma Fourth-Order Derivatives in Mean and Dispersion

Computes the five distinct fourth derivatives of the gamma log-density
with respect to \\\mu\\ and \\\phi\\, in closed form. With \\s =
1/\phi\\ and \\z = y/\mu\\, \$\$\ell^{(\mu\mu\mu\mu)} = \dfrac{s(6 -
24z)}{\mu^4}, \qquad \ell^{(\mu\mu\mu\phi)} = -\dfrac{s^2(6z -
2)}{\mu^3}, \qquad \ell^{(\mu\mu\phi\phi)} = \dfrac{2s^3(1 -
2z)}{\mu^2}, \qquad \ell^{(\mu\phi\phi\phi)} = -\dfrac{6s^4(z -
1)}{\mu},\$\$ and the pure dispersion component follows from the
derivatives in \\s\\ by Faa di Bruno on \\s(\phi) = 1/\phi\\,
\$\$\ell^{(\phi^4)} = f_4 (s')^4 + 6 f_3 (s')^2 s'' + 3 f_2 (s'')^2 + 4
f_2 s' s''' + f_1 s'''',\$\$ with \\s' = -s^2\\, \\s'' = 2s^3\\, \\s'''
= -6s^4\\, \\s'''' = 24s^5\\ and \\f_k\\ the \\k\\th derivative of
\\\ell\\ in \\s\\. The four \\f_k\\ are polygamma functions of \\s\\,
each computed as a polygamma minus its own leading asymptote so that the
digits survive at large shape.

With `expected = TRUE` the expectations are returned, obtained by
replacing \\z\\ with 1 and \\\log z\\ with \\\psi(s) - \log s\\. Both
routes are closed form, so `approx` and `nsim` are ignored.

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

A named list of five numeric vectors, `mu_mu_mu_mu`, `mu_mu_mu_phi`,
`mu_mu_phi_phi`, `mu_phi_phi_phi` and `phi_phi_phi_phi`, each of length
`max(length(y), length(mu), length(phi))`.

## Notation

\\\ell^{(ijkl)}\\ is the fourth derivative of the log-density in
parameters \\i\\, \\j\\, \\k\\ and \\l\\. \\\psi\\ is the digamma
function.

## See also

[`distrib_deriv3.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.Gamma1Distrib.md)
for the order below,
[`distrib_deriv4()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.md)
for the generic, and
[`distrib_expected_hessian.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_expected_hessian.Gamma1Distrib.md)
for the second-order expectation these extend.

## Examples

``` r
d <- gamma1_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, phi = 0.5)
d4 <- distrib_deriv4(d, y, th)
names(d4)
#> [1] "mu_mu_mu_mu"     "mu_mu_mu_phi"    "mu_mu_phi_phi"   "mu_phi_phi_phi" 
#> [5] "phi_phi_phi_phi"

# The four components involving mu, written out.
s <- 1 / 0.5
z <- y / 3
all.equal(d4$mu_mu_mu_mu, s * (6 - 24 * z) / 3^4)
#> [1] TRUE
all.equal(d4$mu_phi_phi_phi, -6 * s^4 * (z - 1) / 3)
#> [1] TRUE

# Expected values: the components odd in z - 1 vanish.
lapply(distrib_deriv4(d, y, th, expected = TRUE), unique)
#> $mu_mu_mu_mu
#> [1] -0.4444444
#> 
#> $mu_mu_mu_phi
#> [1] -0.5925926
#> 
#> $mu_mu_phi_phi
#> [1] -1.777778
#> 
#> $mu_phi_phi_phi
#> [1] 0
#> 
#> $phi_phi_phi_phi
#> [1] -159.6578
#> 

# A central difference of the third order reproduces it.
eps <- 1e-6
up <- distrib_deriv3(d, y, list(mu = 3, phi = 0.5 + eps))$mu_mu_phi
dn <- distrib_deriv3(d, y, list(mu = 3, phi = 0.5 - eps))$mu_mu_phi
all.equal((up - dn) / (2 * eps), d4$mu_mu_phi_phi, tolerance = 1e-4)
#> [1] TRUE
```
