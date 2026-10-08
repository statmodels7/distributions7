# Weibull Moments, Mean and Shape

With \\g_k = \Gamma(1 + k/\sigma)\\, the mean is \\m\\, the variance
\\m^2 (g_2/g_1^2 - 1)\\, the skewness \\(g_3 - 3 g_1 g_2 + 2
g_1^3)/(g_2 - g_1^2)^{3/2}\\ and the excess kurtosis \\(g_4 - 4 g_1
g_3 + 6 g_1^2 g_2 - 3 g_1^4)/(g_2 - g_1^2)^2 - 3\\. The ratios are
formed on the log scale with
[`lgamma()`](https://rdrr.io/r/base/Special.html).

## Arguments

- x:

  A `Weibull3Distrib` object.

- theta:

  A named list with components `mean` and `sigma`.

- ...:

  Unused.

## Value

A numeric vector, one value per parameter row.

## See also

[`weibull3_distrib()`](https://statmodels7.github.io/distributions7/reference/weibull3_distrib.md).

## Examples

``` r
d <- weibull3_distrib()
th <- list(mean = 4, sigma = 1.7)
c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
#> [1] 4.0000000 5.8656899 0.8650234 0.7723790
```
