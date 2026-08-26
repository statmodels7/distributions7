# Excess Kurtosis of the Gamma Distribution in Mean and Dispersion

Closed form: \\\gamma_2 = 6\phi\\, the excess over the Gaussian. It is
linear in the dispersion and free of the mean, and it is exactly
\\3\gamma_1^2/2\\, so a gamma occupies one curve of the
skewness-kurtosis plane whichever parametrization it is written in.

## Arguments

- x:

  A `Gamma1Distrib`, from
  [`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

- theta:

  A named list with components `mu` (positive) and `phi` (positive),
  each a numeric vector of length 1 or `n`. Only `phi` enters the value.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$phi))`, positive throughout.

## Notation

\\\phi \> 0\\ is the dispersion. The mean does not enter a standardized
moment.

## See also

[`skewness.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gamma1Distrib.md),
to which this is tied;
[`kurtosis.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gamma2Distrib.md);
[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md).

## Examples

``` r
d <- gamma1_distrib()

# Six times the dispersion.
all.equal(kurtosis(d, list(mu = 2, phi = 0.25)), 1.5)
#> [1] TRUE

# The family lies on the curve kurtosis = 1.5 skewness^2.
th <- list(mu = 2, phi = 0.25)
all.equal(kurtosis(d, th), 1.5 * skewness(d, th)^2)
#> [1] TRUE
```
