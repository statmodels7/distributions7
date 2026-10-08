# The Scale of a Weibull with a Given Mean

Returns \\b = m / \Gamma(1 + 1/\sigma)\\, formed on the log scale with
[`lgamma()`](https://rdrr.io/r/base/Special.html).

## Usage

``` r
weibull3_scale(theta)
```

## Arguments

- theta:

  A list with the mean and the shape, in that order.

## Value

A numeric vector of scales.

## See also

[`distrib_pdf.Weibull3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Weibull3Distrib.md)
