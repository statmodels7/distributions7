# Excess Kurtosis of the Weibull Distribution

Closed form: \$\$\gamma_2 = \frac{g_4 - 4 g_1 g_3 + 6 g_1^2 g_2 - 3
g_1^4} {(g_2 - g_1^2)^2} - 3, \qquad g_k = \Gamma(1 + k/\sigma).\$\$ The
scale cancels, so the value depends on the shape alone. It is 6 at shape
1, where the family is exponential, dips below zero over a band of
moderate shapes, and climbs again as the shape grows.

## Arguments

- x:

  A `Weibull1Distrib`, from
  [`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

- theta:

  A named list with components `mu` (the scale, positive) and `sigma`
  (the shape, positive), each a numeric vector of length 1 or `n`. The
  fourth moment requires a shape above \\1/4\\.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of excess kurtoses, of length
`max(length(theta$mu), length(theta$sigma))`. Only the shape enters the
value, so a setting that varies the scale alone repeats one number.

## Notation

\\\sigma \> 0\\ is the shape and \\g_k = \Gamma(1 + k/\sigma)\\.

## See also

[`skewness.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Weibull1Distrib.md),
[`variance.Weibull1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Weibull1Distrib.md),
[`weibull_gamma_factors()`](https://statmodels7.github.io/distributions7/reference/weibull_gamma_factors.md),
[`weibull1_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull1_distrib.md).

## Examples

``` r
d <- weibull1_distrib()

# At shape 1 the family is exponential, whose excess kurtosis is 6.
all.equal(kurtosis(d, list(mu = 1, sigma = 1)), 6)
#> [1] TRUE

# It goes negative over a band of moderate shapes.
round(kurtosis(d, list(mu = 1, sigma = c(1, 2, 3.6, 10))), 4)
#> [1]  6.0000  0.2451 -0.2833  0.5702
```
