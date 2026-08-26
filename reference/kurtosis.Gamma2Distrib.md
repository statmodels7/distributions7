# Excess Kurtosis of the Gamma Distribution

Closed form: \\\gamma_2 = 6/a = 6\sigma^2/\mu^2\\ with shape \\a =
\mu^2/\sigma^2\\. It is positive at every parameter value and is exactly
\\3\gamma_1^2/2\\, so the family occupies one curve of the
skewness-kurtosis plane and cannot be tuned along the two axes
separately.

## Arguments

- x:

  A `Gamma2Distrib`, from
  [`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

- theta:

  A named list with components `mu` (positive) and `sigma2` (positive),
  each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$sigma2))`, positive throughout.

## Notation

\\\mu \> 0\\ is the mean, \\\sigma^2 \> 0\\ the variance and \\a =
\mu^2/\sigma^2\\ the shape.

## See also

[`skewness.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gamma2Distrib.md),
to which this is tied;
[`kurtosis.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.ExponentialDistrib.md),
the shape-1 case;
[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

## Examples

``` r
d <- gamma2_distrib()

# The family lies on the curve kurtosis = 1.5 skewness^2.
th <- list(mu = 2, sigma2 = 1)
all.equal(kurtosis(d, th), 1.5 * skewness(d, th)^2)
#> [1] TRUE

# Six over the implied shape.
kurtosis(d, list(mu = 2, sigma2 = c(4, 1, 0.25)))
#> [1] 6.000 1.500 0.375
```
