# Skewness of the Gamma Distribution

Closed form: \\\gamma_1 = 2/\sqrt a = 2\sigma/\mu\\ with shape \\a =
\mu^2/\sigma^2\\. It is positive at every parameter value, the support
being the positive half-line, and it is exactly twice the coefficient of
variation, so a gamma's asymmetry and its relative spread carry the same
information.

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

A numeric vector, of length
`max(length(theta$mu), length(theta$sigma2))`, positive throughout.

## Notation

\\\mu \> 0\\ is the mean, \\\sigma^2 \> 0\\ the variance and \\a =
\mu^2/\sigma^2\\ the shape.

## See also

[`kurtosis.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gamma2Distrib.md),
which is \\3/2\\ times the square of this;
[`variance.Gamma2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gamma2Distrib.md);
[`skewness.ExponentialDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.ExponentialDistrib.md),
the shape-1 case;
[`gamma2_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma2_distrib.md).

## Examples

``` r
d <- gamma2_distrib()

# Twice the coefficient of variation.
all.equal(skewness(d, list(mu = 2, sigma2 = 1)), 2 * sqrt(1) / 2)
#> [1] TRUE

# It flattens as the implied shape grows.
skewness(d, list(mu = 2, sigma2 = c(4, 1, 0.25)))
#> [1] 2.0 1.0 0.5
```
