# The Skewness and Excess Kurtosis of a Weibull Shape

Returns the skewness and the excess kurtosis of a Weibull of shape
\\s\\, which do not depend on the scale, from the ratios \\\Gamma(1 +
k/s)/\Gamma(1 + 1/s)^k\\ formed with
[`lgamma()`](https://rdrr.io/r/base/Special.html).

## Usage

``` r
weibull3_shape_moments(s)
```

## Arguments

- s:

  A numeric vector of shapes.

## Value

A list with components `skew` and `kurt`.

## See also

[`skewness.Weibull3Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Weibull3Distrib.md)
