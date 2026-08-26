# Beta Third-Order Derivatives in Mean and Precision

Computes the four distinct third derivatives of the beta log-density
with respect to \\\mu\\ and \\\phi\\, in closed form as combinations of
\\\psi_2\\, the second derivative of the digamma function, at \\\alpha =
\mu\phi\\, \\\beta = (1-\mu)\phi\\ and \\\phi\\.

**Every component of this order is free of the response.** The data
enter the log-density only through \\(\alpha-1)\log y +
(\beta-1)\log(1-y)\\, which is linear in the shapes and so at most
quadratic in \\(\mu, \phi)\\ through the bilinear map \\\alpha =
\mu\phi\\; three derivatives kill it. The observed and expected values
therefore coincide exactly, and `expected` selects nothing: the same
kernel runs either way and the two results are identical to the bit.
`approx` and `nsim` are ignored for the same reason.

## Arguments

- distrib:

  A `Beta1Distrib` object, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- y:

  A numeric vector of observations. Only its length is used.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `mu` must lie strictly in \\(0, 1)\\ and `phi` must be strictly
  positive.

- expected:

  Logical of length 1, and without effect here, the observed and
  expected third derivatives being the same numbers. Defaults to
  `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). Read by
  the generic, not by this method.

- approx:

  Ignored, a closed form being available.

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
`mu_phi_phi` and `phi_phi_phi`, each of length `length(y)` and constant
within itself when the parameters are. The names enumerate the distinct
multi-indices of order three in two parameters.

## Notation

\\\ell^{(ijk)}\\ is the third derivative of the log-density in
parameters \\i\\, \\j\\ and \\k\\. \\\psi\\ is the digamma function and
\\\psi_m\\ its \\m\\th derivative.

## See also

[`distrib_hessian.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_hessian.Beta1Distrib.md)
for the order below, whose mixed entry is the last one to carry the
response, and
[`distrib_deriv4.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv4.Beta1Distrib.md)
for the order above;
[`distrib_deriv3()`](https://statmodels7.github.io/distributions7/reference/distrib_deriv3.md)
for the generic.

## Examples

``` r
d <- beta1_distrib()
y <- c(0.2, 0.5, 0.8)
th <- list(mu = 0.4, phi = 5)
d3 <- distrib_deriv3(d, y, th)

# Four constants: nothing at this order depends on the observation.
lapply(d3, unique)
#> $mu_mu_mu
#> [1] 31.25
#> 
#> $mu_mu_phi
#> [1] -4.045836
#> 
#> $mu_phi_phi
#> [1] 0.00385982
#> 
#> $phi_phi_phi
#> [1] 0.01036213
#> 

# So asking for the expectation changes nothing.
identical(d3, distrib_deriv3(d, y, th, expected = TRUE))
#> [1] TRUE

# A central difference of the Hessian reproduces the same component.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 0.4 + eps, phi = 5))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, phi = 5))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
#> [1] TRUE
```
