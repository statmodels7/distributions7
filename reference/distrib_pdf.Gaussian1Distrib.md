# Gaussian Probability Density Function

Computes the Gaussian density \$\$f(y; \mu, \sigma) =
\dfrac{1}{\sqrt{2\pi}\\\sigma}
\exp\left\\-\dfrac{1}{2}\left(\dfrac{y-\mu}{\sigma}\right)^2\right\\\$\$
by calling [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html) at
`mean = mu` and `sd = sigma`, so the accuracy and the underflow behavior
are R's own. With `log = TRUE` the logarithm is formed inside
[`dnorm()`](https://rdrr.io/r/stats/Normal.html) and stays finite far
into the tails, where the density itself underflows to zero.

## Arguments

- distrib:

  A `Gaussian1Distrib` object, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

- y:

  A numeric vector of observations. Every real value is in the support,
  so no value is rejected.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive; a zero or negative value
  gives `NaN` with a warning from
  [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html).

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma))`, one value per observation.

## See also

[`distrib_cdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gaussian1Distrib.md)
for the distribution function,
[`distrib_gradient.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian1Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic and
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the family.

## Examples

``` r
d <- gaussian1_distrib()
y <- c(-1.2, 0.3, 2.5)

# The method is stats::dnorm at this parametrization.
all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5)),
          dnorm(y, mean = 0.4, sd = 1.5))
#> [1] TRUE

# A parameter may vary by observation, one value each.
distrib_pdf(d, y, list(mu = c(0, 1, 2), sigma = c(1, 1.5, 2)))
#> [1] 0.1941861 0.2385223 0.1933341

# In the far tail the density underflows and its logarithm does not.
distrib_pdf(d, 40, list(mu = 0, sigma = 1))
#> [1] 0
distrib_pdf(d, 40, list(mu = 0, sigma = 1), log = TRUE)
#> [1] -800.9189
```
