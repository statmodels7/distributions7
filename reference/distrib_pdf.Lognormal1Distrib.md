# Lognormal Probability Density Function

Computes the lognormal density \$\$f(y; \mu, \sigma^2) =
\dfrac{1}{y\sqrt{2\pi\sigma^2}} \exp\left\\-\dfrac{(\log y -
\mu)^2}{2\sigma^2}\right\\, \qquad y \> 0,\$\$ by calling
[`stats::dlnorm()`](https://rdrr.io/r/stats/Lognormal.html) at
`meanlog = mu` and `sdlog = sqrt(sigma2)`. The factor \\1/y\\ is the
Jacobian of the log transformation; it carries no parameter, which is
why every derivative in \\\mu\\ and \\\sigma^2\\ is the Gaussian's read
at \\\log y\\.

With `log = TRUE` the logarithm is formed inside
[`dlnorm()`](https://rdrr.io/r/stats/Lognormal.html) and stays finite
where the density itself underflows.

## Arguments

- distrib:

  A `Lognormal1Distrib` object, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- y:

  A numeric vector of observations. The support is \\(0, \infty)\\; a
  value at or below zero gives 0.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma2` must be strictly positive.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma2))`, one value per observation.

## See also

[`distrib_cdf.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Lognormal1Distrib.md)
for the distribution function,
[`distrib_gradient.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Lognormal1Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gaussian2Distrib.md)
for the law this is the exponential of, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- lognormal1_distrib()
y <- c(0.5, 1.6, 4)
th <- list(mu = 0.5, sigma2 = 0.36)

# The method is stats::dlnorm at meanlog = mu, sdlog = sqrt(sigma2).
all.equal(distrib_pdf(d, y, th), dlnorm(y, 0.5, sqrt(0.36)))
#> [1] TRUE

# Times the Jacobian y, it is the Gaussian density at log y.
all.equal(distrib_pdf(d, y, th) * y,
          distrib_pdf(gaussian2_distrib(), log(y), th))
#> [1] TRUE

# A parameter may vary by observation, one value each.
distrib_pdf(d, y, list(mu = c(0, 0.5, 1), sigma2 = 0.36))
#> [1] 0.6823165 0.4150459 0.1351106

# Far out in the tail the density underflows and its logarithm does not.
distrib_pdf(d, 1e12, th)
#> [1] 0
distrib_pdf(d, 1e12, th, log = TRUE)
#> [1] -1050.39
```
