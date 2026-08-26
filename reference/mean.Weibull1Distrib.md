# Mean of the Weibull Distribution

Closed form: \\E\[Y\] = \mu\\\Gamma(1 + 1/\sigma)\\. The first parameter
of this parametrization is the **scale**, following `gamlss`'s `WEI`, so
it is not the mean; the gamma factor is what separates the two, and it
is 1 only at shape 1, where the family is exponential.

## Arguments

- x:

  A `Weibull1Distrib`, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- theta:

  A named list with components `mu` (the scale, positive) and `sigma`
  (the shape, positive), each a numeric vector of length 1 or `n`. The
  mean diverges as the shape approaches zero.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma))`.

## Details

A mean parametrization was available and was not taken: it would make
every derivative of the family a derivative of the gamma function and of
its inverse. The cost of the choice is paid here, in one gamma
evaluation per setting.

## Notation

\\\mu \> 0\\ is the scale and \\\sigma \> 0\\ the shape, in the
parametrization
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md)
uses.

## See also

[`variance.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Weibull1Distrib.md),
[`skewness.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Weibull1Distrib.md),
[`weibull_gamma_factors()`](https://statmodels7.github.io/distributions7/reference/weibull_gamma_factors.md)
for the shared factors,
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

## Examples

``` r
d <- weibull1_distrib()

# The scale times Gamma(1 + 1 / shape).
all.equal(mean(d, list(mu = 2, sigma = 3)), 2 * gamma(1 + 1 / 3))
#> [1] TRUE

# At shape 1 the family is exponential and the scale is the mean.
mean(d, list(mu = 2, sigma = 1))
#> [1] 2
```
