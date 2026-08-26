# Mean of the Poisson-Inverse Gaussian Distribution

Closed form: \\E\[Y\] = \mu\\, the first cumulant. The family is a
Poisson whose rate is mixed over an inverse Gaussian of mean 1, so the
mixing leaves the mean where the Poisson put it and shows only from the
second cumulant onwards.

## Arguments

- x:

  A `Pig1Distrib`, from
  [`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `sigma`
  (the dispersion, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma))`. The value is the mean
parameter itself, the dispersion entering neither it nor the length.

## Notation

\\\mu \> 0\\ is the mean and \\\sigma \> 0\\ the dispersion of the
mixing inverse Gaussian.

## See also

[`variance.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Pig1Distrib.md),
which exceeds this;
[`mean.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Pig2Distrib.md)
for the orthogonal parametrization;
[`mean.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.PoissonDistrib.md),
the \\\sigma \to 0\\ limit;
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md).

## Examples

``` r
d <- pig1_distrib()

# The first parameter is the mean, and the dispersion does not move it.
mean(d, list(mu = c(1, 2, 3), sigma = 1))
#> [1] 1 2 3
```
