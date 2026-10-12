# The Scale of a Weibull with a Given Mean

`weibull3_scale()` returns \\b = m / \Gamma(1 + 1/\sigma)\\, formed on
the log scale with [`lgamma()`](https://rdrr.io/r/base/Special.html),
and `weibull3_log_scale()` returns \\\log b\\. `weibull3_scale_ok()`
reports where \\b\\ itself is a finite positive number; at a small shape
it underflows to 0, and the density, distribution, quantile and
generator methods then work with \\\log b\\.

## Usage

``` r
weibull3_scale(theta)

weibull3_log_scale(theta)

weibull3_scale_ok(lb)
```

## Arguments

- theta:

  A list with the mean and the shape, in that order.

- lb:

  A numeric vector of log-scales, for `weibull3_scale_ok()`.

## Value

A numeric vector of scales, of log-scales, or of logicals.

## See also

[`distrib_pdf.Weibull3Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Weibull3Distrib.md)
