# Mean of the Skew Normal in the Centered Parametrization

Returns \\\mu\\, the first parameter, which in this parametrization is
the mean by construction. The addition of
[`moment_const()`](https://statmodels7.github.io/distributions7/reference/moment_const.md)
recycles the value to the length the three parameters imply, so a
`theta` whose components vary by observation gets one mean per
observation.

## Arguments

- x:

  A `SkewNormal2Distrib` object, from
  [`skewnormal2_distrib()`](https://statmodels7.github.io/distributions7/reference/skewnormal2_distrib.md).

- theta:

  A named list with components `mu`, `sigma` and `gamma1`, in any order;
  it is aligned here.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of the length the recycled parameters imply.

## See also

[`variance.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.SkewNormal2Distrib.md)
and
[`skewness.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.SkewNormal2Distrib.md),
the other two parameters read back;
[`kurtosis.SkewNormal2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.SkewNormal2Distrib.md),
which is not a parameter; and
[`mean.SkewNormal1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.SkewNormal1Distrib.md)
for the same quantity in the direct parametrization, where it is not.

## Examples

``` r
d <- skewnormal2_distrib()
mean(d, list(mu = 3, sigma = 2, gamma1 = 0.6))
#> [1] 3

# One value per observation when a parameter varies.
mean(d, list(mu = c(0, 3, 7), sigma = 2, gamma1 = 0.6))
#> [1] 0 3 7
```
