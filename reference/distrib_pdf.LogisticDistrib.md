# Logistic Probability Density Function

Computes the logistic density \$\$f(y; \mu, \sigma) =
\dfrac{\exp\left(-\dfrac{y-\mu}{\sigma}\right)}{\sigma \left\[1 +
\exp\left(-\dfrac{y-\mu}{\sigma}\right)\right\]^2}\$\$ by calling
[`stats::dlogis()`](https://rdrr.io/r/stats/Logistic.html) at
`location = mu` and `scale = sigma`. The density is symmetric about
\\\mu\\ and its tails decay exponentially, like
\\e^{-\|y-\mu\|/\sigma}\\, so they are heavier than a Gaussian's.

## Arguments

- distrib:

  A `LogisticDistrib` object, from
  [`logistic_distrib()`](https://statmodels7.github.io/distributions7/reference/logistic_distrib.md).

- y:

  A numeric vector of observations. Every real value is in the support,
  so no value is rejected.

- theta:

  A named list with components `mu` and `sigma`, each a numeric vector
  of length 1 or of the length of `y`. A component of length 1 is
  recycled. `sigma` must be strictly positive; a zero or negative value
  gives `NaN` with a warning from
  [`stats::dlogis()`](https://rdrr.io/r/stats/Logistic.html).

- log:

  Logical of length 1. When `TRUE` the log-density is returned. Defaults
  to `FALSE`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of densities, of length
`max(length(y), length(mu), length(sigma))`, one value per observation.

## See also

[`distrib_cdf.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_cdf.LogisticDistrib.md)
for the distribution function,
[`distrib_gradient.LogisticDistrib()`](https://statmodels7.github.io/distributions7/reference/distrib_gradient.LogisticDistrib.md)
for the derivatives of the log-density, and
[`distrib_pdf()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.md)
for the generic.

## Examples

``` r
d <- logistic_distrib()
y <- c(-1.2, 0.3, 2.5)

# The method is stats::dlogis at this parametrization.
all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5)),
          dlogis(y, location = 0.4, scale = 1.5))
#> [1] TRUE

# Symmetric about mu: equal densities at equal distances either side.
distrib_pdf(d, 0.4 + c(-2, 2), list(mu = 0.4, sigma = 1.5))
#> [1] 0.1100607 0.1100607

# The density and the distribution function are linked by
# f = F (1 - F) / sigma, the logistic sigmoid's own derivative.
F <- distrib_cdf(d, y, list(mu = 0.4, sigma = 1.5))
all.equal(distrib_pdf(d, y, list(mu = 0.4, sigma = 1.5)), F * (1 - F) / 1.5)
#> [1] TRUE
```
