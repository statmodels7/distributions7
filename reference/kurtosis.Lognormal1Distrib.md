# Excess Kurtosis of the Lognormal Distribution

Closed form: \$\$\gamma_2 = e^{4\sigma^2} + 2e^{3\sigma^2} +
3e^{2\sigma^2} - 6,\$\$ free of \\\mu\\. The leading term grows like
\\e^{4\sigma^2}\\, so the tail weight explodes with the spread of the
logarithm: at \\\sigma^2 = 1\\ it is 110.9, and the numerical route
would need a quadrature far into the tail to see it at all.

## Arguments

- x:

  A `Lognormal1Distrib`, from
  [`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

- theta:

  A named list with components `mu` (any real value) and `sigma2`
  (positive), each a numeric vector of length 1 or `n`. Only `sigma2`
  enters the value; the leading exponential overflows to `Inf` past
  \\\sigma^2 \approx 177\\.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$sigma2))`, positive throughout.

## Notation

\\\sigma^2 \> 0\\ is the variance of \\\log Y\\.

## See also

[`skewness.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Lognormal1Distrib.md),
[`variance.Lognormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Lognormal1Distrib.md),
[`lognormal1_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal1_distrib.md).

## Examples

``` r
d <- lognormal1_distrib()

# The published form, written out.
all.equal(kurtosis(d, list(mu = 0, sigma2 = 1)),
          exp(4) + 2 * exp(3) + 3 * exp(2) - 6)
#> [1] TRUE

# It climbs steeply with the spread of the logarithm.
round(kurtosis(d, list(mu = 0, sigma2 = c(0.1, 0.5, 1))), 3)
#> [1]   1.856  18.507 110.936
```
