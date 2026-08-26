# Mean of the Negative Binomial Distribution

Closed form, replacing the numerical default: \\E\[Y\] = \mu\\. This
parametrization carries the mean as its first parameter, so the method
reads it off and recycles it to the length the two parameters imply.
Nothing is integrated and nothing is summed.

## Arguments

- x:

  A `NegBin2Distrib`, from
  [`negbin2_distrib()`](https://statmodels7.github.io/distributions7/reference/negbin2_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `theta`
  (the dispersion, positive), each a numeric vector of length 1 or `n`.
  Aligned and validated by name, so a missing or out-of-bounds component
  throws.

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
