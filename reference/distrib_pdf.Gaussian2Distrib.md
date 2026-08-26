# Gaussian Probability Density Function in Mean and Variance

Computes the Gaussian density \$\$f(y; \mu, \sigma^2) =
\dfrac{1}{\sqrt{2\pi\sigma^2}}
\exp\left\\-\dfrac{(y-\mu)^2}{2\sigma^2}\right\\\$\$ by calling
[`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html) at `mean = mu`
and `sd = sqrt(sigma2)`, so the accuracy and the underflow behavior are
R's own. With `log = TRUE` the logarithm is formed inside
[`dnorm()`](https://rdrr.io/r/stats/Normal.html) and stays finite far
into the tails, where the density itself underflows to zero.

## Arguments

- distrib:

  A `Gaussian2Distrib` object, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- y:

  A numeric vector of observations. Every real value is in the support,
  so no value is rejected.

- theta:

  A named list with components `mu` and `sigma2`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma2` must be strictly positive; a zero or negative value
  makes `sqrt(sigma2)` zero or `NaN` and
  [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html) warns.

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma2))`, one value per observation.

## See also

[`distrib_cdf.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gaussian2Distrib.md)
for the distribution function,
[`distrib_gradient.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian2Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gaussian1Distrib.md)
for the same density in the standard deviation, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- gaussian2_distrib()
y <- c(-1.2, 0.3, 2.5)

# The method is stats::dnorm at the square root of the variance.
all.equal(distrib_pdf(d, y, list(mu = 1, sigma2 = 4)),
          dnorm(y, mean = 1, sd = 2))
#> [1] TRUE

# Same law as gaussian1 at sigma = sqrt(sigma2), to the last bit.
all.equal(distrib_pdf(d, y, list(mu = 1, sigma2 = 4)),
          distrib_pdf(gaussian1_distrib(), y, list(mu = 1, sigma = 2)))
#> [1] TRUE

# A parameter may vary by observation, one value each.
distrib_pdf(d, y, list(mu = c(0, 1, 2), sigma2 = c(1, 2.25, 4)))
#> [1] 0.1941861 0.2385223 0.1933341

# In the far tail the density underflows and its logarithm does not.
distrib_pdf(d, 40, list(mu = 0, sigma2 = 1))
#> [1] 0
distrib_pdf(d, 40, list(mu = 0, sigma2 = 1), log = TRUE)
#> [1] -800.9189
```
