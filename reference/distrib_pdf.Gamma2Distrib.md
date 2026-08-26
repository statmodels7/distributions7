# Gamma Probability Density Function in Mean and Variance

Computes the gamma density \$\$f(y; \mu, \sigma^2) =
\dfrac{\lambda^{\alpha}}{\Gamma(\alpha)} y^{\alpha - 1} e^{-\lambda y},
\qquad \alpha = \dfrac{\mu^2}{\sigma^2}, \quad \lambda =
\dfrac{\mu}{\sigma^2}, \quad y \> 0,\$\$ by calling
[`stats::dgamma()`](https://rdrr.io/r/stats/GammaDist.html) at that
shape and rate. With `log = TRUE` the logarithm is formed inside
[`dgamma()`](https://rdrr.io/r/stats/GammaDist.html) and stays finite
where the density itself underflows.

The density is unbounded at the origin when \\\sigma^2 \> \mu^2\\, where
the shape falls below one; it is flat there when \\\sigma^2 = \mu^2\\,
the exponential case, and vanishes at the origin when \\\sigma^2 \<
\mu^2\\.

## Arguments

- distrib:

  A `Gamma2Distrib` object, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- y:

  A numeric vector of observations. The support is \\(0, \infty)\\; a
  negative value gives 0.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. Both must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma2))`, one value per observation.

## See also

[`distrib_cdf.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gamma2Distrib.md)
for the distribution function,
[`distrib_gradient.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gamma2Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gamma1Distrib.md)
for the same density in the dispersion, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- gamma2_distrib()
y <- c(1, 3, 5)
th <- list(mu = 3, sigma2 = 2)

# The method is stats::dgamma at shape mu^2/sigma2 and rate mu/sigma2.
all.equal(distrib_pdf(d, y, th),
          dgamma(y, shape = 9 / 2, rate = 3 / 2))
#> [1] TRUE

# The same law as gamma1 at phi = sigma2/mu^2.
all.equal(distrib_pdf(d, y, th),
          distrib_pdf(gamma1_distrib(), y, list(mu = 3, phi = 2 / 9)))
#> [1] TRUE

# A parameter may vary by observation, one value each.
distrib_pdf(d, y, list(mu = c(1, 3, 9), sigma2 = 2))
#> [1] 0.241970725 0.276927214 0.001525943

# Far out in the tail the density underflows and its logarithm does not.
distrib_pdf(d, 1e4, th)
#> [1] 0
distrib_pdf(d, 1e4, th, log = TRUE)
#> [1] -14968.39
```
