# Mean of the Gaussian Distribution

Closed form: \\E\[Y\] = \mu\\. The first parameter of this
parametrization is the mean itself, so the method reads it off and
recycles it to the length the two parameters imply.

## Arguments

- x:

  A `Gaussian1Distrib`, from
  [`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md).

- theta:

  A named list with components `mu` (the mean, any real value) and
  `sigma` (the standard deviation, positive), each a numeric vector of
  length 1 or `n`. Aligned and validated by name, so a missing or
  out-of-bounds component throws.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma))`.

## See also

[`variance.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gaussian1Distrib.md),
[`skewness.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.Gaussian1Distrib.md)
and
[`kurtosis.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.Gaussian1Distrib.md),
which are the family's other three moments;
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
for the family.

## Examples

``` r
d <- gaussian1_distrib()

# The first parameter is the mean, and the scale does not move it.
mean(d, list(mu = c(-1, 0, 4), sigma = 3))
#> [1] -1  0  4
```
