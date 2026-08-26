# The Cauchy Distribution Has No Mean

Returns `NaN`. The Cauchy density decays like \\y^{-2}\\, so \\\int
\|y\| f(y)\\\mathrm{d}y\\ diverges and the first moment does not exist.
`NaN` is returned directly, so that no quadrature is attempted: a
numerical integration over a divergent integral returns whatever its
truncation gives, a number that moves with the panel layout and reads
like an estimate.

## Arguments

- x:

  A `CauchyDistrib`, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- theta:

  A named list with components `mu` and `sigma` (positive), each a
  numeric vector of length 1 or `n`. The values are not read, only their
  lengths, and the list is still aligned and validated, so a missing or
  out-of-bounds component throws as it would for a family whose moments
  exist.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `NaN`, of length
`max(length(theta$mu), length(theta$sigma))`.

## Details

What the family does have is a **median**, equal to \\\mu\\, and a
half-interquartile range, equal to \\\sigma\\. Both are exact at every
parameter value and both come from
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md).
A Cauchy sample mean is itself Cauchy with the same scale, so it does
not settle as the sample grows and cannot stand in for the missing
moment either.

## Notation

\\\mu\\ is the location and \\\sigma \> 0\\ the scale, in the
parametrization
[`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md)
uses. Neither is a moment.

## See also

[`skewness.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/skewness.CauchyDistrib.md)
for the full argument;
[`variance.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.CauchyDistrib.md)
and
[`kurtosis.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.CauchyDistrib.md),
`NaN` for the same reason;
[`distrib_quantile()`](https://statmodels7.github.io/distributions7/reference/distrib_quantile.md)
for the median and the quartiles.

## Examples

``` r
d <- cauchy_distrib()

# No moment exists, so this is NaN at every parameter value.
mean(d, list(mu = c(0, 1, 2), sigma = 1))
#> [1] NaN NaN NaN

# The median is mu and the half-interquartile range is sigma.
distrib_quantile(d, c(0.25, 0.5, 0.75), list(mu = 3, sigma = 2))
#> [1] 1 3 5
```
