# The Cauchy Distribution Has No Kurtosis

Returns `NaN`. The excess kurtosis is the fourth standardized central
moment less three, so it needs the first four moments of \\Y\\, and the
Cauchy has none of them: its density decays like \\y^{-2}\\, and \\\int
\|y\|^p f(y)\\\mathrm{d}y\\ diverges already at \\p = 1\\. `NaN` is
returned directly, so that the method does not attempt a quadrature.

## Arguments

- x:

  A `CauchyDistrib`, from
  [`cauchy_distrib()`](https://statmodels7.github.io/distributions7/reference/cauchy_distrib.md).

- theta:

  A named list with components `mu` and `sigma` (positive), in which
  each component is a numeric vector of length 1 or `n`. The values are
  not read, only their lengths, and the list is still aligned and
  validated.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of `NaN`, of length
`max(length(theta$mu), length(theta$sigma))`.

## Details

A large positive excess kurtosis usually indicates heavy tails. The
Cauchy has heavy tails, but its fourth moment does not exist, so the
excess kurtosis cannot measure their weight. A Student t with \\\nu \>
4\\ is the nearest family that reports a number, and \\6/(\nu-4)\\ grows
without bound as \\\nu\\ falls towards 4.

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
[`variance.CauchyDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.CauchyDistrib.md),
`NaN` for the same reason;
[`kurtosis.StudentT1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.StudentT1Distrib.md),
which reports a number above four degrees of freedom.

## Examples

``` r
d <- cauchy_distrib()

# The value is NaN at every parameter value.
kurtosis(d, list(mu = 0, sigma = c(1, 2)))
#> [1] NaN NaN

# The Student t is the family that does report a number, above nu = 4.
kurtosis(student_t1_distrib(), list(mu = 0, sigma = 1, nu = c(1, 5, 10)))
#> [1] NaN   6   1
```
