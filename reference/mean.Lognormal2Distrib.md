# Lognormal Moments, Mean and Variance of Y

The mean is \\m\\ and the variance \\v\\. With \\e = v/m^2\\ the
skewness is \\(e + 3)\sqrt{e}\\ and the excess kurtosis \\e(16 + 15e +
6e^2 + e^3)\\, which are \\(w + 2)\sqrt{w - 1}\\ and \\w^4 + 2w^3 +
3w^2 - 6\\ at \\w = 1 + e\\ written without the subtraction that loses
digits when \\e\\ is small.

## Arguments

- x:

  A `Lognormal2Distrib` object.

- theta:

  A named list with components `mean` and `var`.

- ...:

  Unused.

## Value

A numeric vector, one value per parameter row.

## See also

[`lognormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/lognormal2_distrib.md).

## Examples

``` r
d <- lognormal2_distrib()
th <- list(mean = 3, var = 2)
c(mean(d, th), variance(d, th), skewness(d, th), kurtosis(d, th))
#> [1] 3.000000 2.000000 1.518970 4.364579
```
