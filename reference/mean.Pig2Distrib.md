# Mean of the Orthogonal Poisson-Inverse Gaussian Distribution

Closed form: \\E\[Y\] = \mu\\. This is the same law as
[`pig1_distrib()`](https://statmodels7.github.io/distributions7/reference/pig1_distrib.md)
written by a dispersion \\\alpha\\ chosen so that the mean and the
dispersion are orthogonal in the expected information; the mapping is
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md),
and the mean is the first parameter either way.

## Arguments

- x:

  A `Pig2Distrib`, from
  [`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

- theta:

  A named list with components `mu` (the mean, positive) and `alpha`
  (the orthogonal dispersion, positive), each a numeric vector of length
  1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$alpha))`. The value is the mean
parameter itself, the dispersion entering neither it nor the length.

## Notation

\\\mu \> 0\\ is the mean and \\\alpha \> 0\\ the orthogonal dispersion.

## See also

[`variance.Pig2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Pig2Distrib.md),
where the dispersion does enter;
[`pig2_sigma()`](https://statmodels7.github.io/distributions7/reference/pig2_sigma.md)
for the mapping;
[`mean.Pig1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Pig1Distrib.md);
[`pig2_distrib()`](https://statmodels7.github.io/distributions7/reference/pig2_distrib.md).

## Examples

``` r
d <- pig2_distrib()

# The first parameter is the mean, and the dispersion does not move it.
mean(d, list(mu = c(1, 2, 3), alpha = 1))
#> [1] 1 2 3
```
