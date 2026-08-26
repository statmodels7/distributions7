# Skewness of the Gamma Distribution in Mean and Dispersion

Closed form: \\\gamma_1 = 2\sqrt\phi\\. The mean cancels, so the
asymmetry is carried by the dispersion alone: it is twice the
coefficient of variation and vanishes as the dispersion goes to zero,
where the family approaches a Gaussian.

## Arguments

- x:

  A `Gamma1Distrib`, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

- theta:

  A named list with components `mu` (positive) and `phi` (positive),
  each a numeric vector of length 1 or `n`. Only `phi` enters the value,
  and the result is recycled to the length both imply.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector, of length `max(length(theta$mu), length(theta$phi))`,
positive throughout.

## Notation

\\\phi \> 0\\ is the dispersion. The mean does not enter a standardized
moment.

## See also

[`kurtosis.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gamma1Distrib.md),
which is \\3/2\\ times the square of this;
[`skewness.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gamma2Distrib.md),
the same quantity in the mean-variance parametrization;
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

## Examples

``` r
d <- gamma1_distrib()

# Twice the square root of the dispersion, whatever the mean.
skewness(d, list(mu = c(1, 2, 3), phi = 0.25))
#> [1] 1 1 1

# It agrees with the mean-variance parametrization on the same law.
all.equal(skewness(d, list(mu = 2, phi = 0.25)),
          skewness(gamma2_distrib(), list(mu = 2, sigma2 = 1)))
#> [1] TRUE
```
