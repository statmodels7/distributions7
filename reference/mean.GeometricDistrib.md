# Mean of the Geometric Distribution

Closed form: \\E\[Y\] = \mu\\. The family is parametrized by its mean,
the expected number of failures before the first success, and the
success probability it implies is \\1/(1+\mu)\\.

## Arguments

- x:

  A `GeometricDistrib`, from
  [`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

- theta:

  A named list with one component, `mu` (the mean, positive), a numeric
  vector of length 1 or `n`.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, the length of `theta$mu`.

## See also

[`variance.GeometricDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.GeometricDistrib.md),
which exceeds this;
[`mean.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.NegBin2Distrib.md),
of which this is the \\\theta = 1\\ case;
[`geometric_distrib()`](https://statmodels7.github.io/distributions7/reference/geometric_distrib.md).

## Examples

``` r
d <- geometric_distrib()

# The one parameter is the mean.
mean(d, list(mu = c(0.5, 3, 20)))
#> [1]  0.5  3.0 20.0
```
