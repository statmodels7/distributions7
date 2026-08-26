# Beta Probability Density Function in the Shapes

Computes the beta density \$\$f(y; \alpha, \beta) =
\dfrac{y^{\alpha-1}(1-y)^{\beta-1}} {B(\alpha, \beta)}, \qquad 0 \< y \<
1,\$\$ with \\B\\ the beta function, by calling
[`stats::dbeta()`](https://rdrr.io/r/stats/Beta.html) at
`shape1 = alpha` and `shape2 = beta`. With `log = TRUE` the logarithm is
formed inside [`dbeta()`](https://rdrr.io/r/stats/Beta.html) and stays
finite where the density itself underflows.

The density is unbounded at 0 when \\\alpha \< 1\\ and at 1 when \\\beta
\< 1\\; at \\\alpha = \beta = 1\\ it is the uniform.

## Arguments

- distrib:

  A `Beta2Distrib` object, from
  [`beta2_distrib()`](https://statmodels7.github.io/distributions7/reference/beta2_distrib.md).

- y:

  A numeric vector of observations in \\(0, 1)\\. A value outside \\\[0,
  1\]\\ gives 0, and the endpoints give 0, a finite value or `Inf`
  according to the shapes.

- theta:

  A named list with components `alpha` and `beta`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(alpha), length(beta))`, one value per
observation.

## See also

[`distrib_cdf.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Beta2Distrib.md)
for the distribution function,
[`distrib_gradient.Beta2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Beta2Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf.Beta1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Beta1Distrib.md)
for the same density in the mean and the precision, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- beta2_distrib()
y <- c(0.1, 0.3, 0.7)
th <- list(alpha = 2, beta = 5)

# The method is stats::dbeta at these two shapes.
all.equal(distrib_pdf(d, y, th), dbeta(y, 2, 5))
#> [1] TRUE

# The same law as beta1 at mu = alpha/(alpha + beta), phi = alpha + beta.
all.equal(distrib_pdf(d, y, th),
          distrib_pdf(beta1_distrib(), y, list(mu = 2 / 7, phi = 7)))
#> [1] TRUE

# Both shapes 1 is the uniform.
distrib_pdf(d, y, list(alpha = 1, beta = 1))
#> [1] 1 1 1

# Near the boundary the density underflows and its logarithm does not.
distrib_pdf(d, 1e-40, list(alpha = 20, beta = 5))
#> [1] 0
distrib_pdf(d, 1e-40, list(alpha = 20, beta = 5), log = TRUE)
#> [1] -1737.698
```
