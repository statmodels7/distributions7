# Gaussian Probability Density Function in Mean and Precision

Computes the Gaussian density \$\$f(y; \mu, \tau) =
\sqrt{\dfrac{\tau}{2\pi}}
\exp\left\\-\dfrac{\tau(y-\mu)^2}{2}\right\\\$\$ by calling
[`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html) at `mean = mu`
and `sd = 1/sqrt(tau)`, so the accuracy and the underflow behavior are
R's own. With `log = TRUE` the logarithm is formed inside
[`dnorm()`](https://rdrr.io/r/stats/Normal.html) and stays finite far
into the tails, where the density itself underflows to zero.

## Arguments

- distrib:

  A `Gaussian3Distrib` object, from
  [`gaussian3_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian3_distrib.md).

- y:

  A numeric vector of observations. Every real value is in the support,
  so no value is rejected.

- theta:

  A named list with components `mu` and `tau`, each a numeric vector of
  length 1 or of the length of `y`. A component of length 1 is recycled.
  `tau` must be strictly positive; a zero gives `sd = Inf` and a density
  of 0, and a negative value gives `NaN` with a warning from
  [`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html).

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(tau))`, one value per observation.

## See also

[`distrib_cdf.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.Gaussian3Distrib.md)
for the distribution function,
[`distrib_gradient.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.Gaussian3Distrib.md)
for the derivatives of the log-density,
[`distrib_pdf.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gaussian1Distrib.md)
for the same density in the standard deviation, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- gaussian3_distrib()
y <- c(-1.2, 0.3, 2.5)

# The method is stats::dnorm at the reciprocal square root of tau.
all.equal(distrib_pdf(d, y, list(mu = 1, tau = 0.25)),
          dnorm(y, mean = 1, sd = 2))
#> [1] TRUE

# Same law as gaussian1 at sigma = 1/sqrt(tau), to the last bit.
all.equal(distrib_pdf(d, y, list(mu = 1, tau = 0.25)),
          distrib_pdf(gaussian1_distrib(), y, list(mu = 1, sigma = 2)))
#> [1] TRUE

# A parameter may vary by observation, one value each.
distrib_pdf(d, y, list(mu = c(0, 1, 2), tau = c(1, 0.25, 0.0625)))
#> [1] 0.19418605 0.18762017 0.09895942

# In the far tail the density underflows and its logarithm does not.
distrib_pdf(d, 40, list(mu = 0, tau = 1))
#> [1] 0
distrib_pdf(d, 40, list(mu = 0, tau = 1), log = TRUE)
#> [1] -800.9189
```
