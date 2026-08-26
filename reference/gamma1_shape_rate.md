# The Shape and Rate a Mean and Dispersion Imply

Converts the `gamma1` parameters to the shape and rate
[`stats::dgamma()`](https://rdrr.io/r/stats/GammaDist.html) takes: \\a =
1/\phi\\ and \\b = 1/(\phi\mu)\\. Under that pairing the mean \\a/b\\ is
\\\mu\\ and the variance \\a/b^2\\ is \\\phi\mu^2\\. Nothing is
validated; a non-positive `phi` or `mu` propagates as an infinite or
negative value to the caller.

## Usage

``` r
gamma1_shape_rate(theta)
```

## Arguments

- theta:

  A named list with components `mu` and `phi`, each a numeric vector.
  The two are used elementwise, so components of different lengths
  recycle in the usual way.

## Value

A list of two numeric vectors, `shape` and `rate`, of the lengths the
arithmetic produces.

## See also

[`gamma1_distrib()`](https://statmodels7.github.io/distributions7/reference/gamma1_distrib.md)
for the family and
[`distrib_pdf.Gamma1Distrib()`](https://statmodels7.github.io/distributions7/reference/distrib_pdf.Gamma1Distrib.md)
for the density this feeds.
