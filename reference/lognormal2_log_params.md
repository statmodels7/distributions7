# The Log-Scale Parameters of a Mean and a Variance

Returns the mean and the standard deviation of \\\log Y\\ for the mean
\\m\\ and the variance \\v\\ of \\Y\\: \\S = \log(1 + v/m^2)\\, formed
by [`log1p()`](https://rdrr.io/r/base/Log.html), `sdlog` \\\sqrt{S}\\
and `meanlog` \\\log m - S/2\\.

## Usage

``` r
lognormal2_log_params(theta)
```

## Arguments

- theta:

  A list with the mean and the variance, in that order.

## Value

A list with components `meanlog` and `sdlog`.

## See also

[`distrib_pdf.Lognormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Lognormal2Distrib.md)
