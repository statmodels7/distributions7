# Mean of the Negative Binomial Distribution

The closed form that replaces the numerical default is \\E\[Y\] = \mu\\.
This parametrization carries the mean as its first parameter, so the
method reads it off and recycles it to the length that the two
parameters imply. The method performs no integration and no summation.

## Arguments

- x:

  A `NegBin2Distrib`, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `theta`
  (the dispersion, positive); each is a numeric vector of length 1 or
  `n`. It is aligned and validated by name, so a missing or
  out-of-bounds component signals an error.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$theta))`.

## See also

[`variance.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.NegBin2Distrib.md),
which exceeds this by \\\mu^2/\theta\\;
[`skewness.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/skewness.NegBin2Distrib.md)
and
[`kurtosis.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/kurtosis.NegBin2Distrib.md);
[`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md)
for the family.

## Examples

``` r
d <- negbin2_distrib()

# The mean is the first parameter, and the dispersion does not move it.
mean(d, list(mu = c(1, 2, 3), theta = 2))
#> [1] 1 2 3
```
