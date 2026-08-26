# Mean of the Gaussian Distribution in Location and Variance

Closed form: \\E\[Y\] = \mu\\. This is the same law as
[`gaussian1_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian1_distrib.md)
written by its variance, and the mean is the location in all three of
the toolkit's Gaussian parametrizations.

## Arguments

- x:

  A `Gaussian2Distrib`, from
  [`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

- theta:

  A named list with components `mu` (the mean, any real value) and
  `sigma2` (the variance, positive), each a numeric vector of length 1
  or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$sigma2))`.

## See also

[`variance.Gaussian2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Gaussian2Distrib.md),
where the parametrizations differ;
[`mean.Gaussian1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gaussian1Distrib.md)
and
[`mean.Gaussian3Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.Gaussian3Distrib.md);
[`gaussian2_distrib()`](https://statmodels7.github.io/distributions7/reference/gaussian2_distrib.md).

## Examples

``` r
d <- gaussian2_distrib()

# The first parameter is the mean, and the variance does not move it.
mean(d, list(mu = c(-1, 0, 4), sigma2 = 9))
#> [1] -1  0  4
```
