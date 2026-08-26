# Mean of the Poisson Distribution

Closed form: \\E\[Y\] = \mu\\. The family carries one parameter and it
is the mean, so this method reads it off. It is also the variance, which
is the equidispersion the family is defined by.

## Arguments

- x:

  A `PoissonDistrib`, from
  [`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

- theta:

  A named list with one component, `mu` (the mean, positive), a numeric
  vector of length 1 or `n`. Aligned and validated by name, so a
  non-positive value throws.

- ...:

  Unused, and accepted so that the signature matches the generic's.

## Value

A numeric vector of means, the length of `theta$mu`.

## See also

[`variance.PoissonDistrib()`](https://statmodels7.github.io/distributions7/reference/variance.PoissonDistrib.md),
equal to this;
[`mean.NegBin2Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.NegBin2Distrib.md)
and
[`mean.NegBin1Distrib()`](https://statmodels7.github.io/distributions7/reference/mean.NegBin1Distrib.md),
the overdispersed count families;
[`poisson_distrib()`](https://statmodels7.github.io/distributions7/reference/poisson_distrib.md).

## Examples

``` r
d <- poisson_distrib()

# The one parameter is the mean, and it is also the variance.
c(mean(d, list(mu = 3)), variance(d, list(mu = 3)))
#> [1] 3 3
```
