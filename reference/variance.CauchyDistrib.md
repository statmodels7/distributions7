# The Cauchy Distribution Has No Variance

Returns `NaN`. The variance is the second central moment, and it needs a
first moment to center on; the Cauchy has neither. Its density decays
like \\y^{-2}\\, so \\\int \|y\|^p f(y)\\\mathrm{d}y\\ diverges for
every \\p \ge 1\\, and the failure at \\p = 2\\ is the worse of the two.
`NaN` is returned directly, so that no quadrature is attempted.

## Arguments

- x:

  A `CauchyDistrib`, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- theta:

  A named list with components `mu` and `sigma` (positive), each a
  numeric vector of length 1 or `n`. The values are not read, only their
  lengths, and the list is still aligned and validated.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `NaN`, of length
`max(length(theta$mu), length(theta$sigma))`.
[`std_dev()`](https://statmodels7.github.io/distributions7/reference/std_dev.md)
gives `NaN` too, being the square root of this.

## Details

The spread a Cauchy does have is its half-interquartile range, exactly
\\\sigma\\, available from
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md).
Note that \\\sigma\\ is a legitimate scale parameter: it sets the width
of the density in the ordinary way, and only its identification with a
standard deviation fails.

## Notation

\\\mu\\ is the location and \\\sigma \> 0\\ the scale, in the
parametrization
[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
uses.

## See also

[`skewness.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.CauchyDistrib.md)
for the full argument;
[`mean.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.CauchyDistrib.md)
and
[`kurtosis.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.CauchyDistrib.md),
`NaN` for the same reason;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the quartiles;
[`variance.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.StudentT1Distrib.md),
which is this family at \\\nu = 1\\ and reports `NaN` there too.

## Examples

``` r
d <- cauchy_distrib()

# NaN at every parameter value, and so is the standard deviation.
c(variance(d, list(mu = 0, sigma = 2)), std_dev(d, list(mu = 0, sigma = 2)))
#> [1] NaN NaN

# The half-interquartile range is sigma, and it is exact.
q <- distrib_quantile(d, c(0.25, 0.75), list(mu = 3, sigma = 2))
diff(q) / 2
#> [1] 2
```
