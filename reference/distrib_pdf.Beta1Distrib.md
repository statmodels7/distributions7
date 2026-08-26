# Beta Probability Density Function in Mean and Precision

Computes the beta density \$\$f(y; \mu, \phi) =
\dfrac{\Gamma(\phi)}{\Gamma(\mu\phi)\\ \Gamma((1-\mu)\phi)}\\
y^{\mu\phi - 1}(1-y)^{(1-\mu)\phi - 1}, \qquad 0 \< y \< 1,\$\$ by
calling [`stats::dbeta()`](https://rdrr.io/r/stats/Beta.html) at
`shape1 = mu * phi` and `shape2 = (1 - mu) * phi`. With `log = TRUE` the
logarithm is formed inside
[`dbeta()`](https://rdrr.io/r/stats/Beta.html) and stays finite where
the density itself underflows.

The shape of the density is governed by whether each shape parameter
exceeds one. It is unbounded at 0 when \\\mu\phi \< 1\\ and at 1 when
\\(1-\mu)\phi \< 1\\; at \\\mu = 1/2\\ and \\\phi = 2\\ both shapes are
1 and the density is the uniform.

## Arguments

- distrib:

  A `Beta1Distrib` object, from
  [`beta1_distrib()`](https://statmodels7.github.io/distributions7/reference/beta1_distrib.md).

- y:

  A numeric vector of observations. The support is \\(0, 1)\\; a value
  outside \\\[0, 1\]\\ gives 0, and the endpoints give 0, a finite value
  or `Inf` according to the shapes.

- theta:

  A named list with components `mu` and `phi`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `mu` must lie strictly in \\(0, 1)\\ and `phi` must be strictly
  positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(phi))`, one value per observation.

## See also

[`distrib_cdf.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Beta1Distrib.md)
for the distribution function,
[`distrib_gradient.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Beta1Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Beta2Distrib.md)
for the same density in the variance, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- beta1_distrib()
y <- c(0.2, 0.5, 0.8)
th <- list(mu = 0.4, phi = 5)

# The method is stats::dbeta at shape1 = mu phi and shape2 = (1 - mu) phi.
all.equal(distrib_pdf(d, y, th), dbeta(y, 0.4 * 5, 0.6 * 5))
#> [1] TRUE

# Both shapes are 1 at mu = 1/2, phi = 2, where the beta is the uniform.
distrib_pdf(d, y, list(mu = 0.5, phi = 2))
#> [1] 1 1 1

# A parameter may vary by observation, one value each.
distrib_pdf(d, y, list(mu = c(0.2, 0.5, 0.8), phi = 5))
#> [1] 2.048000 1.697653 2.048000

# Near the boundary the density underflows and its logarithm does not.
distrib_pdf(d, 1e-40, list(mu = 0.4, phi = 50))
#> [1] 0
distrib_pdf(d, 1e-40, list(mu = 0.4, phi = 50), log = TRUE)
#> [1] -1715.996
```
