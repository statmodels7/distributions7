# Mean of the Laplace Distribution in Location and Rate

Closed form, replacing the numerical default: \\E\[Y\] = \mu\\. This is
the same family as
[`laplace_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace_distrib.md)
written by its rate \\\lambda = 1/\sigma\\, and the mean does not see
the difference: it is the location in either parametrization.

## Arguments

- x:

  A `Laplace2Distrib`, from
  [`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

- theta:

  A named list with components `mu` (the location) and `lambda` (the
  rate, positive), each a numeric vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, of length
`max(length(theta$mu), length(theta$lambda))`.

## See also

[`variance.Laplace2Distrib()`](https://statmodels7.github.io/distributions7/reference/variance.Laplace2Distrib.md),
where the parametrization does show;
[`mean.LaplaceDistrib()`](https://statmodels7.github.io/distributions7/reference/mean.LaplaceDistrib.md)
for the scale form;
[`laplace2_distrib()`](https://statmodels7.github.io/distributions7/reference/laplace2_distrib.md).

## Examples

``` r
d <- laplace2_distrib()

# The location is the mean, and the rate does not move it.
mean(d, list(mu = c(-1, 0, 4), lambda = 2))
#> [1] -1  0  4
```
