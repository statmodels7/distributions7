# Beta Third-Order Derivatives in Mean and Precision

Computes the four distinct third derivatives of the beta log-density
with respect to \\\mu\\ and \\\phi\\, in closed form as combinations of
\\\psi_1\\ and \\\psi_2\\, the first and second derivatives of the
digamma function, at \\\alpha = \mu\phi\\, \\\beta = (1-\mu)\phi\\ and
\\\phi\\.

Every component of this order is free of the response. The data enter
the log-density only through \\(\alpha-1)\log y + (\beta-1)\log(1-y)\\,
which is linear in the shapes and so at most quadratic in \\(\mu,
\phi)\\ through the bilinear map \\\alpha = \mu\phi\\, and its third
derivatives are zero. On the parameter scale the observed and expected
values therefore coincide, and the method runs the same kernel for
either value of `expected`. On the link scale the generic adds terms in
the score and the Hessian, which depend on the response, so the two
differ there.

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

  Logical of length 1, without effect here because the observed and
  expected third derivatives are the same numbers. Defaults to `FALSE`.

- scale:

  One of `"parameter"` (the default) or `"link"`, matched by
  [`base::match.arg()`](https://rdrr.io/r/base/match.arg.html). The
  generic reads it, and this method does not.

- method:

  `NULL` or an object inheriting from
  [ExpectMethod](https://statmodels7.github.io/distributions7/reference/ExpectMethod.md).
  The generic applies an object before dispatching on the family, so
  this method receives `NULL`, and it returns the closed form of the
  family. See
  [expected_derivative_methods](https://statmodels7.github.io/distributions7/reference/expected_derivative_methods.md).

- ...:

  Unused, and accepted so that the signature matches the generic's.

- threads:

  A single positive integer, how many threads the kernel may use. Below
  a fixed number of observations, set per kernel, the kernel runs
  sequentially. Defaults to `1L`.

## Value

A named list of four numeric vectors, `mu_mu_mu`, `mu_mu_phi`,
`mu_phi_phi` and `phi_phi_phi`, each of length
`max(length(y), length(mu), length(phi))` and constant within itself
when the parameters are. The names enumerate the distinct multi-indices
of order three in two parameters.

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

# Four constants: every component at this order is free of the observation.
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

# So on the parameter scale expected = TRUE changes nothing.
identical(d3, distrib_deriv3(d, y, th, expected = TRUE))
#> [1] TRUE

# A central difference of the Hessian reproduces the same component.
eps <- 1e-6
up <- distrib_hessian(d, y, list(mu = 0.4 + eps, phi = 5))$mu_mu
dn <- distrib_hessian(d, y, list(mu = 0.4 - eps, phi = 5))$mu_mu
all.equal((up - dn) / (2 * eps), d3$mu_mu_mu, tolerance = 1e-4)
#> [1] TRUE
```
