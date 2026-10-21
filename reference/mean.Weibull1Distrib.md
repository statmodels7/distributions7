# Mean of the Weibull Distribution

The closed form is \\E\[Y\] = \mu\\\Gamma(1 + 1/\sigma)\\. The first
parameter of this parametrization is the scale, following `gamlss`'s
`WEI`, so it is not the mean. The two differ by the gamma factor, which
equals 1 only at shape 1, where the family is exponential.

## Arguments

- x:

  A `Weibull1Distrib`, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- theta:

  A named list with components `mu` (the scale, positive) and `sigma`
  (the shape, positive); each is a numeric vector of length 1 or `n`.
  The mean diverges as the shape approaches zero.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma))`.

## Details

The package also ships the Weibull in its mean and shape,
[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md).
In the scale parametrization of this page the mean needs one gamma
evaluation per setting.

## Notation

\\\mu \> 0\\ is the scale and \\\sigma \> 0\\ the shape, in the
parametrization that
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
uses.

## See also

[`variance.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Weibull1Distrib.md),
[`skewness.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Weibull1Distrib.md),
[`weibull_moment_terms()`](https://statmodels7.github.io/distributions7/reference/weibull_moment_terms.md)
for the terms in the shape,
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

## Examples

``` r
d <- weibull1_distrib()

# The mean is the scale times Gamma(1 + 1 / shape).
all.equal(mean(d, list(mu = 2, sigma = 3)), 2 * gamma(1 + 1 / 3))
#> [1] TRUE

# At shape 1 the family is exponential and the scale is the mean.
mean(d, list(mu = 2, sigma = 1))
#> [1] 2
```
